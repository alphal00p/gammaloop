use std::collections::BTreeMap;

use color_eyre::Result;
use eyre::{WrapErr, eyre};
use gammaloop_tracing_filter::debug_instrument;
use idenso::{
    color::{ColorSimplifier, ColorSimplifySettings},
    dirac::{AGS, GammaSimplifier, GammaSimplifySettings},
    representations::Bispinor,
    shorthands::{
        UndoShorthands,
        chain::Chain,
        metric::MetricSimplifier,
        schoonschip::{Schoonschip, SchoonschipSettings},
    },
    tensor::{SymbolicNetExt, SymbolicNetParse},
};

use linnet::half_edge::{
    HedgeGraph,
    builder::HedgeGraphBuilder,
    involution::{EdgeIndex, HedgePair},
    subgraph::{ModifySubSet, SuBitGraph, SubGraphLike, SubSetLike, SubSetOps},
};
use spenso::{
    network::{
        graph::NetworkEdge,
        library::symbolic::ETS,
        parsing::{
            ParseSettings, SchoonschipExpansionMode, ShorthandParsing, StructureInferenceMode,
        },
        tags::SPENSO_TAG,
    },
    shadowing::TensorCollectExt,
    structure::{
        representation::{Minkowski, RepName},
        slot::{DummyAind, IsAbstractSlot, ParseableAind, Slot},
    },
};
use symbolica::{
    atom::{AliasedAtom, Atom, AtomCore, AtomView, FunctionBuilder},
    domains::{atom::AtomField, integer::Z, rational::Q},
    function,
    id::Replacement,
    parse, parse_lit,
    poly::{PolyVariable, series::Series},
    symbol,
};
use symbolica_utils::ReplaceBuilderExt;
use vakint::{Vakint, VakintExpression, vakint_symbol};

use crate::{
    debug_tags,
    graph::{LMBext, LoopMomentumBasis},
    numerator::aind::Aind,
    utils::{GS, W_},
    uv::{
        ApproximationType, UltravioletGraph,
        approx::{
            ForestNodeLike, Rooted, UVCtx,
            local_4d::{FourDSector, Local4dCts},
        },
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
        // The enclosing forest Taylor operation rescales the vacuum mass in
        // this coefficient along with its active momenta. The bookkeeping
        // factor restores only consumed loop measures; physical projections
        // and the next integration set it to one.
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
    if integrand.contains_symbol(AGS.sigma) {
        return Err(eyre!(
            "Sigma tensors are not supported by d-dimensional analytic UV numerator algebra"
        ));
    }
    // Only the analytically integrated UV subgraph reaches this function.
    // Expose spin/vector slots before closing traces; the representation-specific
    // collectors below retain factorized scalar coefficients.
    let integrand =
        if integrand.contains_symbol(AGS.gamma) || integrand.contains_symbol(SPENSO_TAG.trace) {
            // Scalar products already bind their contracted slots. Materializing
            // a powered dot before symbolic multiplication reuses its dummy index
            // and turns (p.k)^2 into p^2 k^2. Protect completed scalar spectators,
            // while keeping spin/vector contractions with open slots visible.
            // The symbolic network also checks hidden shorthand
            // slots without requiring concrete tensor dimensions or execution.
            let mut used = integrand.get_all_symbols(true);
            let mut aliases = BTreeMap::<Atom, Atom>::new();
            let mut serial = 0usize;
            let protected = integrand.replace_map(|arg, _context, out| {
                // The complete scalar expression in a production denominator is
                // metadata for propagator conversion, not a tensor-network factor.
                // Protect it alongside scalar spectators before closing spin traces.
                if matches!(arg, AtomView::Fun(fun) if fun.get_symbol() == GS.den)
                    || (matches!(arg, AtomView::Fun(fun)
                    if fun.get_symbol() == SPENSO_TAG.dot || fun.get_symbol() == ETS.metric)
                        && !arg.contains_symbol(AGS.gamma)
                        && !arg.contains_symbol(SPENSO_TAG.trace)
                        && arg
                            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                                shorthand_parsing: ShorthandParsing::expand_all(),
                                ..Default::default()
                            })
                            .is_ok_and(|network| network.graph.dangling_indices().is_empty()))
                {
                    let alias = aliases.entry(arg.to_owned()).or_insert_with(|| {
                        loop {
                            let symbol = symbol!(&format!("gammaloop::uv_scalar_product_{serial}"));
                            serial += 1;
                            if used.insert(symbol) {
                                break Atom::var(symbol);
                            }
                        }
                    });
                    **out = alias.clone();
                }
            });
            let mut protected = AliasedAtom::from(protected);
            for (original, alias) in aliases {
                protected.register_alias(alias, original);
            }
            let prepared = protected
                .get_root()
                .parse_to_symbolic_net::<Aind>(&ParseSettings {
                    shorthand_parsing: ShorthandParsing::Expand {
                        schoonschip: SchoonschipExpansionMode::full(),
                        trace: false,
                        chain: true,
                    },
                    ..Default::default()
                })
                .map_err(|error| eyre!("invalid analytic UV spin tensor notation: {error}"))?
                .simple_execute::<()>()?;
            protected.map_root(|_| prepared).into_inner()
        } else {
            integrand.clone()
        };
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
    // Keep each Taylor term's propagator powers. Vakint receives the individual
    // denominator topologies: collecting factors across the sum clears denominators
    // and manufactures higher-rank numerator factors before reduction.
    // Color is a spectator of the spin algebra; collecting all chains and traces
    // here distributes its shared coefficient across the kinematic terms.
    let collected = schoonschip.simplify_metrics().collect_gamma_chains();
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

    if dotted
        .replace(function!(
            SPENSO_TAG.trace,
            Bispinor {}.to_symbolic([W_.d_]),
            W_.x___
        ))
        .matches()
    {
        return Err(eyre!(
            "an unresolved Dirac trace remains after d-dimensional UV numerator algebra"
        ));
    }
    if dotted
        .replace(function!(
            SPENSO_TAG.trace,
            Minkowski {}.to_symbolic([W_.d_]),
            W_.x___
        ))
        .matches()
    {
        return Err(eyre!(
            "an unresolved Lorentz trace remains after d-dimensional UV numerator algebra"
        ));
    }
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
            ApproximationType::MUV | ApproximationType::PolePart | ApproximationType::IR => {
                // Integrate the complete signed operator in each retained frame.
                // Linearity permits this sector sum; no soft branch is discarded
                // before the nested forest has supplied its completed coefficient.
                let mut integrated = Atom::Zero;
                for sector in integrand.active_sectors() {
                    let physical = sector
                        .atom
                        .replace(GS.integrated_loop_scale)
                        .with(Atom::one());
                    let simplified = simplify(&physical)?;
                    if !simplified.is_zero() {
                        integrated += self
                            .integrate(&simplified, sector, ctx, current, given)
                            .wrap_err_with(|| {
                                format!(
                                    "while integrating {scheme} counterterm {} with retained component frames {:?}",
                                    current.subgraph().string_label(),
                                    sector.active_components,
                                )
                            })?;
                    }
                }
                let marker = UvMarker::new(ctx.settings);
                let integrated = marker.apply(
                    UvOperation::Integrate,
                    marker_current.subgraph(),
                    marker_given.subgraph(),
                    &integrated,
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
        sector: &FourDSector,
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
        let integrand_vakint = to_vakint_integrand(
            integrand,
            graph,
            current.subgraph(),
            given.subgraph(),
            &sector.active_components,
            &settings.vakint,
            false,
        )?;

        let propagator = function!(
            vakint::symbols::S.prop,
            W_.a_,
            W_.b_,
            W_.mom_,
            W_.mass_,
            W_.e_
        )
        .to_pattern();
        let retained_loop_count = sector
            .active_components
            .iter()
            .flat_map(|(_, _, lmb)| lmb.loop_edges.iter())
            .collect::<std::collections::HashSet<_>>()
            .len();
        let mut res = Atom::Zero;
        for source_term in integrand_vakint.0 {
            // A product of independent one-loop denominators can carry a
            // coupled numerator. Integrate one factor at a time, treating all
            // spectator loops as fixed vectors for its angular projection.
            let mut loop_factors = std::collections::BTreeMap::<i64, Atom>::new();
            let mut masses = std::collections::HashSet::new();
            let mut independent = true;
            for matched in source_term.integral.pattern_match(&propagator, None, None) {
                let momentum = &matched[&W_.mom_];
                let loop_id = [1, -1].into_iter().find_map(|sign| {
                    let momentum = momentum * Atom::num(sign);
                    let AtomView::Fun(vector) = momentum.as_view() else {
                        return None;
                    };
                    (vector.get_symbol() == vakint::symbols::S.k && vector.get_nargs() == 1)
                        .then(|| i64::try_from(vector.get(0)).ok())
                        .flatten()
                });
                let Some(loop_id) = loop_id else {
                    independent = false;
                    break;
                };
                masses.insert(matched[&W_.mass_].clone());
                *loop_factors.entry(loop_id).or_insert_with(Atom::one) *= propagator
                    .replace_wildcards(&matched)
                    .map_err(|error| eyre!(error))?;
            }
            let factorized = independent
                && retained_loop_count > 1
                && loop_factors.len() == retained_loop_count
                && masses.len() > 1;
            let loop_ids = loop_factors.keys().copied().collect::<Vec<_>>();
            let stages = if factorized {
                loop_factors
                    .into_iter()
                    .map(|(loop_id, propagators)| {
                        (
                            function!(vakint::symbols::S.topo, propagators),
                            loop_ids
                                .iter()
                                .copied()
                                .filter(|id| *id != loop_id)
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect::<Vec<_>>()
            } else {
                vec![(source_term.integral.clone(), Vec::new())]
            };
            let mut coefficient = source_term.numerator.clone();
            if factorized {
                // Preserve even custom normalization hooks that are not a
                // product of identical per-loop factors. The converter's
                // additional normalization already belongs to the full measure.
                let normalization = self
                    .vakint_settings
                    .get_integral_normalization_factor_atom()?;
                coefficient *= normalization
                    .replace(vakint::symbols::S.n_loops.to_pattern())
                    .with(Atom::num(retained_loop_count))
                    / normalization
                        .replace(vakint::symbols::S.n_loops.to_pattern())
                        .with(Atom::one())
                        .pow(retained_loop_count);
            }
            // Keep the same conservative epsilon target at every stage. Later
            // vacuum poles can consume positive orders, while Vakint reserves
            // coefficient and normalization poles independently.
            for (integral, spectators) in stages {
                let stage_loop_count = if factorized { 1 } else { retained_loop_count };
                let external = function!(vakint::symbols::S.p, W_.i_, W_.x___).to_pattern();
                let mut used_ids = coefficient
                    .pattern_match(&external, None, None)
                    .map(|matched| i64::try_from(matched[&W_.i_].as_view()))
                    .collect::<Result<std::collections::HashSet<_>, _>>()
                    .map_err(|error| eyre!(error))?;
                let mut next_id = 0i64;
                let mut restore_spectators = Vec::new();
                for spectator in spectators {
                    while used_ids.contains(&next_id) {
                        next_id = next_id.checked_add(1).ok_or_else(|| {
                            eyre!("no free external momentum ID for vacuum spectator")
                        })?;
                    }
                    used_ids.insert(next_id);
                    let source = function!(vakint::symbols::S.k, spectator, W_.x___);
                    let target = function!(vakint::symbols::S.p, next_id, W_.x___);
                    coefficient = coefficient
                        .replace(source.to_pattern())
                        .with(target.to_pattern());
                    restore_spectators.push(Replacement::new(target, source));
                }
                let mut stage_term = source_term.clone();
                stage_term.integral = integral;
                stage_term.numerator = coefficient;
                let mut integrand_vakint = VakintExpression(vec![stage_term]);
                let mut same_momentum = std::collections::BTreeMap::<Atom, Vec<_>>::new();
                for matched in integrand_vakint.0[0]
                    .integral
                    .pattern_match(&propagator, None, None)
                {
                    let momentum = &matched[&W_.mom_];
                    same_momentum
                        .entry(momentum.clone().min(-momentum))
                        .or_default()
                        .push((
                            usize::try_from(matched[&W_.a_].as_view())
                                .map_err(|error| eyre!(error))?,
                            matched[&W_.mass_].clone(),
                            matched[&W_.e_].clone(),
                        ));
                }
                for (momentum, propagators) in same_momentum {
                    let (_, first_mass, _) = &propagators[0];
                    if propagators.iter().all(|(_, mass, _)| mass == first_mass)
                        || momentum.contains_symbol(vakint::symbols::S.p)
                        || propagators.iter().any(|(_, mass, power)| {
                            mass.contains_symbol(W_.x_)
                                || !i64::try_from(power.as_view()).is_ok_and(|power| power > 0)
                        })
                    {
                        continue;
                    }

                    // A mixed-mass one-loop vacuum is a sum of single-mass
                    // tadpoles: 1/(D_a D_b) = (1/D_a - 1/D_b)/(m_a²-m_b²).
                    // The same identity applies to any equal-momentum group
                    // with other vacuum propagators held fixed. Partial-fraction
                    // only the denominator, including raised powers; the graph
                    // numerator remains an unchanged factor.
                    let loop_square = Atom::var(W_.x_);
                    let denominator = propagators
                        .iter()
                        .fold(Atom::one(), |product, (_, mass, power)| {
                            product * (&loop_square - mass).pow(-power)
                        });
                    let fractions = denominator
                        .try_to_rational_polynomial::<_, _, u16>(&Q, &Z, [W_.x_])?
                        .apart_factored_denominators(0);
                    for term in std::mem::take(&mut integrand_vakint.0) {
                        let powers = term
                            .integral
                            .pattern_match(&propagator, None, None)
                            .map(|matched| {
                                Ok((
                                    usize::try_from(matched[&W_.a_].as_view())
                                        .map_err(|error| eyre!(error))?,
                                    matched[&W_.e_].clone(),
                                ))
                            })
                            .collect::<Result<std::collections::BTreeMap<_, _>>>()?;
                        let graph = vakint::graph::Graph::new_from_atom(
                            term.integral.as_view(),
                            *powers.keys().next_back().unwrap(),
                        )?;
                        for (coefficient, denominator, power) in &fractions {
                            let denominator = denominator.to_expression();
                            let constant = denominator.replace(W_.x_).with(Atom::Zero);
                            let leading = denominator.replace(W_.x_).with(Atom::one()) - &constant;
                            let (survivor, _, _) = propagators.iter().find(|(_, mass, _)| {
                                (&denominator - &leading * (&loop_square - mass)).expand().is_zero()
                            }).ok_or_else(|| eyre!("vacuum partial fraction has an unexpected denominator {denominator}"))?;
                            let coefficient = coefficient.to_expression() / leading.pow(*power);
                            if coefficient.contains_symbol(W_.x_) {
                                return Err(eyre!(
                                    "vacuum partial fraction retained loop dependence in its coefficient"
                                ));
                            }
                            let mut contracted = graph.clone();
                            for (id, _, _) in &propagators {
                                if id != survivor {
                                    contracted.contract_edges(&[*id].into_iter().collect());
                                }
                            }
                            // Removing a serial line contracts its endpoints.
                            // Dropping it without that contraction would silently
                            // lower a sunset's graph rank although both loop
                            // generators still occur in its denominators.
                            if contracted.to_symbolica_graph(false).num_loops() != stage_loop_count
                            {
                                return Err(eyre!(
                                    "vacuum partial fraction changed the retained {stage_loop_count}-loop integration domain for component {reduced_label}"
                                ));
                            }
                            let mut fraction = term.clone();
                            fraction.numerator *= coefficient;
                            fraction.integral = function!(
                                vakint::symbols::S.topo,
                                contracted
                                    .edges
                                    .values()
                                    .fold(Atom::one(), |product, edge| {
                                        product
                                            * function!(
                                                vakint::symbols::S.prop,
                                                edge.id,
                                                function!(
                                                    vakint::symbols::S.edge,
                                                    edge.left_node_id,
                                                    edge.right_node_id
                                                ),
                                                &edge.momentum,
                                                &edge.mass,
                                                if edge.id == *survivor {
                                                    Atom::num(*power)
                                                } else {
                                                    powers[&edge.id].clone()
                                                }
                                            )
                                    })
                            );
                            integrand_vakint.0.push(fraction);
                        }
                    }
                }
                // Only a vacuum topology whose actual masses all vanish is scaleless.
                // Certify polynomial loop dependence as well: an external scale hidden
                // inside a nonpolynomial numerator must not be erased with the soft jet.
                integrand_vakint.0.retain(|term| {
                    let propagators = term
                        .integral
                        .pattern_match(&propagator, None, None)
                        .map(|matched| {
                            matched[&W_.mass_].is_zero()
                                && matched[&W_.mom_].contains_symbol(vakint::symbols::S.k)
                        })
                        .collect::<Vec<_>>();
                    if propagators.is_empty()
                        || propagators.iter().any(|massless_loop| !massless_loop)
                        || term.integral.contains_symbol(vakint::symbols::S.p)
                    {
                        return true;
                    }
                    let mut pending = vec![term.numerator.as_view()];
                    let mut scaleless = true;
                    while let Some(view) = pending.pop() {
                        if !view.contains_symbol(vakint::symbols::S.k) {
                            continue;
                        }
                        match view {
                            AtomView::Add(sum) => pending.extend(sum.iter()),
                            AtomView::Mul(product) => pending.extend(product.iter()),
                            AtomView::Pow(power) => {
                                let (base, exponent) = power.get_base_exp();
                                if i64::try_from(exponent).is_ok_and(|power| power >= 0) {
                                    pending.push(base);
                                } else {
                                    scaleless = false;
                                    break;
                                }
                            }
                            AtomView::Fun(function)
                                if function.get_symbol() == vakint::symbols::S.k
                                    && function.get_nargs() > 0
                                    && usize::try_from(function.get(0)).is_ok() => {}
                            AtomView::Fun(function)
                                if function.get_symbol() == vakint::symbols::S.dot
                                    && function.get_nargs() == 2
                                    && function.iter().all(|argument| {
                                        matches!(argument, AtomView::Fun(vector)
                                            if [vakint::symbols::S.k, vakint::symbols::S.p]
                                                .contains(&vector.get_symbol())
                                                && vector.get_nargs() == 1
                                                && usize::try_from(vector.get(0)).is_ok())
                                    }) => {}
                            _ => {
                                scaleless = false;
                                break;
                            }
                        }
                    }
                    if scaleless {
                        debug_tags!(#uv, #integrated, #vakint, #trace;
                            stage = "certified_scaleless_vacuum",
                            log.integral = term.integral,
                            log.numerator = term.numerator,
                            "Dropping a massless vacuum term with polynomial loop numerator"
                        );
                    }
                    !scaleless
                });
                if integrand_vakint.0.is_empty() {
                    coefficient = Atom::Zero;
                    break;
                }

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

                // The analytic backends identify a vacuum mass through a
                // symbol or its square. Alias other exact nonzero mass values
                // only in propagator mass slots, then restore the values below.
                let mut mass_aliases = std::collections::BTreeMap::new();
                let mut used_symbols = integrand_vakint
                    .0
                    .iter()
                    .flat_map(|term| [&term.integral, &term.numerator])
                    .flat_map(|atom| atom.get_all_symbols(true))
                    .collect::<std::collections::HashSet<_>>();
                let mut mass_serial = 0;
                for term in &mut integrand_vakint.0 {
                    let mut replacements = Vec::new();
                    for matched in term.integral.pattern_match(&propagator, None, None) {
                        let mass = &matched[&W_.mass_];
                        if mass.contains_symbol(vakint::symbols::S.k)
                            || mass.contains_symbol(vakint::symbols::S.p)
                        {
                            return Err(eyre!(
                                "a vacuum propagator mass must be independent of momentum"
                            ));
                        }
                        if mass.is_zero()
                            || matches!(mass.as_view(), AtomView::Var(_))
                            || matches!(mass.as_view(), AtomView::Pow(power)
                                if matches!(power.get_base_exp().0, AtomView::Var(_))
                                    && power.get_base_exp().1 == Atom::num(2).as_view())
                        {
                            continue;
                        }
                        let alias = mass_aliases.entry(mass.clone()).or_insert_with(|| {
                            loop {
                                let symbol = vakint_symbol!(format!(
                                    "integrated_mass_squared_{mass_serial}"
                                ));
                                mass_serial += 1;
                                if used_symbols.insert(symbol) {
                                    break Atom::var(symbol);
                                }
                            }
                        });
                        replacements.push(Replacement::new(
                            propagator
                                .replace_wildcards(&matched)
                                .map_err(|error| eyre!(error))?,
                            function!(
                                vakint::symbols::S.prop,
                                &matched[&W_.a_],
                                &matched[&W_.b_],
                                &matched[&W_.mom_],
                                &*alias,
                                &matched[&W_.e_]
                            ),
                        ));
                    }
                    term.integral = term.integral.replace_multiple(&replacements);
                }
                integrand_vakint.canonicalize(
                    self.vakint_settings,
                    &self.vakint.topologies,
                    false,
                )?;
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
                for (term_index, term) in integrand_vakint.0.iter_mut().enumerate() {
                    term.numerator =
                        self.simplify_projected_numerator(&term.numerator, current.topo_order())?;
                    debug_tags!(#uv, #integrated, #vakint, #trace;
                        stage = "projected_numerator_after_d_dimensional_algebra",
                        term_index = %term_index,
                        log.numerator = term.numerator,
                        "Completed projected numerator algebra before Laurent expansion"
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

                coefficient = Atom::from(integrand_vakint).replace_multiple(
                    mass_aliases
                        .into_iter()
                        .map(|(mass, alias)| Replacement::new(alias, mass))
                        .collect::<Vec<_>>(),
                );
                if coefficient.contains_symbol(vakint::symbols::S.k) {
                    return Err(eyre!(
                        "integrated vacuum factor retained an active loop momentum"
                    ));
                }
                coefficient = coefficient.replace_multiple(&restore_spectators);
            }
            res += coefficient;
        }

        debug_tags!(#uv,#integrated,#vakint,#trace,#raw;
            log.res = res,
            "Raw post vakint "
        );

        let mut res = Self::restore_numerator(res, current.topo_order());

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

        // println!("\nIntegrated CT:\n{}\n", res);
        // Normalize the fully integrated UV subgraph before Laurent projection
        // and reinsertion: identical scalar masters can arrive as 2*(A+B) or
        // 2*A+2*B. No cograph numerator is present at this analytic boundary.
        Ok(res.expand())
    }

    fn simplify_projected_numerator(
        &self,
        numerator: &Atom,
        topology_order: usize,
    ) -> Result<Atom> {
        // Tensor projection can identify Lorentz slots belonging to different
        // gamma matrices. Complete their d-dimensional algebra before Vakint
        // expands scalar coefficients in epsilon; otherwise an evanescent
        // contraction can multiply a pole only after its finite term was lost.
        let numerator = Vakint::convert_to_dot_notation(self.vakint_settings, numerator.as_view())?;
        let numerator = Self::restore_numerator(numerator, topology_order);
        let numerator = simplify(&numerator)?;
        Self::ensure_resolved_lorentz_contractions(&numerator)?;
        let numerator = Self::dimensionally_regularized(&numerator)
            .undo_schoonschip::<Aind>()?
            .undo_chain::<Aind>()?
            .undo_trace::<Aind>()?
            .metric_shorthand_to_dot();
        Ok(Self::to_vakint_numerator(&numerator))
    }

    fn ensure_resolved_lorentz_contractions(numerator: &Atom) -> Result<()> {
        let minkowski = Minkowski {}.to_symbolic([GS.dim]).get_symbol().unwrap();
        if !numerator.contains_symbol(minkowski) {
            return Ok(());
        }
        // Inspect each expanded analytic branch using the tensor parser's
        // incidence rules. Expose compact inter-chain contractions first;
        // keeping only non-momentum Lorentz tensors allows scalar products and
        // slashed external momenta while detecting unresolved tensor operators.
        let explicit = numerator
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::full(),
                    trace: false,
                    chain: true,
                },
                ..Default::default()
            })
            .map_err(|error| eyre!("invalid analytic UV Lorentz tensor notation: {error}"))?
            .simple_execute::<()>()?
            .expand();
        let terms = if let AtomView::Add(sum) = explicit.as_view() {
            sum.iter().collect::<Vec<_>>()
        } else {
            vec![explicit.as_view()]
        };
        let settings = ParseSettings {
            shorthand_parsing: ShorthandParsing::Opaque {
                inference: StructureInferenceMode::Fast,
            },
            ..Default::default()
        };
        for term in terms {
            let mut unresolved_inner_product = false;
            let tensors = term.replace_map(|view, _, out| {
                if !matches!(view, AtomView::Num(_)) && !view.contains_symbol(minkowski) {
                    // Remove scalar subtrees atomically: replacing the d in
                    // 1/(d-1) alone would manufacture a singular coefficient.
                    // Numerical exponents of tensor powers keep their value.
                    **out = Atom::one();
                    return;
                }
                match view {
                    AtomView::Fun(f) => {
                        if [SPENSO_TAG.dot, ETS.metric].contains(&f.get_symbol())
                            && view.contains_symbol(minkowski)
                            && !f
                                .iter()
                                .all(|arg| Slot::<Minkowski, Aind>::try_from(arg).is_ok())
                        {
                            // Full shorthand parsing must expose every compact
                            // Lorentz contraction before scalar coefficients can
                            // be discarded by this diagnostic.
                            unresolved_inner_product = true;
                        }
                        let momentum = [
                            GS.emr_mom,
                            GS.loop_mom,
                            vakint::symbols::S.p,
                            vakint::symbols::S.k,
                        ]
                        .contains(&f.get_symbol());
                        **out = if !momentum
                            && f.iter()
                                .any(|arg| Slot::<Minkowski, Aind>::try_from(arg).is_ok())
                        {
                            view.to_owned()
                        } else {
                            Atom::one()
                        };
                    }
                    AtomView::Var(_) => **out = Atom::one(),
                    _ => {}
                }
            });
            if unresolved_inner_product {
                return Err(eyre!(
                    "an unresolved compact Lorentz inner product remains after d-dimensional UV numerator algebra"
                ));
            }
            let network = tensors
                .parse_to_symbolic_net::<Aind>(&settings)
                .map_err(|error| {
                    eyre!("invalid residual d-dimensional tensor structure: {error}")
                })?;
            let external = network.graph.dangling_indices();
            for (pair, _, edge) in network.graph.graph.iter_edges() {
                if pair.is_paired()
                    && let NetworkEdge::Slot(slot) = edge.data
                    && slot.rep_name().symbol() == minkowski
                    && !external.contains(slot)
                {
                    return Err(eyre!(
                        "an unresolved internal Lorentz tensor contraction remains after d-dimensional UV numerator algebra; a d-dimensional tensor-operator reduction is required before analytic integration"
                    ));
                }
            }
        }
        Ok(())
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

    fn restore_numerator(mut res: Atom, topology_order: usize) -> Atom {
        // Projector indices may close two opaque tensor leaves. Restore their
        // complete Minkowski representation before exposing the tensor bodies.
        let mink = Minkowski {}.new_rep(GS.dim);
        let mut projected_indices = std::collections::BTreeMap::new();
        res = res.replace_map(|view, _, out| {
            if matches!(view, AtomView::Fun(index) if index.get_symbol() == vakint::symbols::S.tensor_index) {
                let replacement = projected_indices.entry(view.to_owned()).or_insert_with(|| loop {
                    let candidate = mink.to_symbolic([Aind::new_dummy().to_atom()]);
                    if res.pattern_match(&candidate.to_pattern(), None, None).next().is_none() {
                        break candidate;
                    }
                });
                **out = replacement.clone();
            }
        });

        res = res
            // Opaque tensor slots are needed only while Vakint validates and
            // projects the Lorentz domain. Restore the caller's tensor leaves.
            .replace(function!(vakint::symbols::S.tensor, W_.x_, W_.x___))
            .with(Atom::var(W_.x_))
            .replace(parse_lit!(vakint::cl2))
            .with(parse_lit!(cl2))
            .replace(parse_lit!(vakint::sqrt3))
            .with(parse_lit!(sqrt(3)));

        let vk_metric = vakint_symbol!("g");

        // apply metric
        res = res
            .replace(vakint::symbols::S.p.call_args([W_.i_, W_.j_]))
            .when(W_.j_.filter(|r| r.is_integer().is_true()))
            .with(
                vakint::symbols::S.p.call_args([
                    Atom::var(W_.i_),
                    mink.to_symbolic([GS
                        .uvaind
                        .call_args([Atom::num(topology_order), Atom::var(W_.j_)])]),
                ]),
            )
            .replace(
                vakint::symbols::S
                    .p
                    .call_args([Atom::var(W_.i_), vakint::symbols::S.dot_dummy_ind(W_.j_)]),
            )
            .when(W_.j_.filter(|r| r.is_integer().is_true()))
            .with(
                vakint::symbols::S.p.call_args([
                    Atom::var(W_.i_),
                    mink.to_symbolic([GS
                        .uvaind
                        .call_args([Atom::num(topology_order), Atom::var(W_.j_)])]),
                ]),
            )
            .replace(vakint::symbols::S.p.call_args([W_.x__]))
            .with(GS.emr_mom.call_args([W_.x__]));
        res = res
            .replace(function!(vk_metric, W_.x_, W_.y_) * function!(GS.emr_mom, W_.x___, W_.x_))
            .with(function!(GS.emr_mom, W_.x___, W_.y_))
            .replace(function!(
                vk_metric,
                vakint::symbols::S.dot_dummy_ind(W_.x_),
                W_.y_
            ))
            .when(W_.x_.filter(|r| r.is_integer().is_true()))
            .with(function!(
                vk_metric,
                mink.to_symbolic([GS
                    .uvaind
                    .call_args([Atom::num(topology_order), Atom::var(W_.x_)])]),
                W_.y_
            ))
            .replace(function!(
                vk_metric,
                W_.x_,
                vakint::symbols::S.dot_dummy_ind(W_.y_)
            ))
            .when(W_.y_.filter(|r| r.is_integer().is_true()))
            .with(function!(
                vk_metric,
                mink.to_symbolic([GS
                    .uvaind
                    .call_args([Atom::num(topology_order), Atom::var(W_.y_)])]),
                W_.x_
            ))
            .replace(function!(vk_metric, W_.x_, W_.y_))
            .with(function!(ETS.metric, W_.x_, W_.y_));

        res = res.replace(vakint::symbols::S.cmplx_i).with(Atom::i());

        res
    }

    fn to_vakint_numerator(numerator: &Atom) -> Atom {
        let numerator = numerator
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
        // Preserve the Lorentz slots of arbitrary spin tensors without
        // teaching Vakint Spenso's representation syntax or tensor algebra.
        // Direct slots use the same parser as the tensor network. The original
        // leaf stays opaque, including its spin and color indices.
        numerator.replace_map(|view, _, out| {
            let AtomView::Fun(tensor) = view else { return };
            if [
                vakint::symbols::S.k,
                vakint::symbols::S.p,
                vakint::symbols::S.g,
                vakint::symbols::S.dot,
                vakint::symbols::S.tensor,
            ]
            .contains(&tensor.get_symbol())
            {
                **out = view.to_owned();
                return;
            }
            let slots = tensor
                .iter()
                .filter(|argument| Slot::<Minkowski, Aind>::try_from(*argument).is_ok())
                .collect::<Vec<_>>();
            if !slots.is_empty() {
                let mut wrapped = FunctionBuilder::new(vakint::symbols::S.tensor).add_arg(view);
                for slot in slots {
                    wrapped = wrapped.add_arg(slot);
                }
                **out = wrapped.finish();
            }
        })
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
    active_components: &[(SuBitGraph, SuBitGraph, LoopMomentumBasis)],
    settings: &VakintSettings,
    substitute_masses_to_m_uv: bool,
) -> Result<VakintExpression> {
    let reduced_label = reduced.string_label();
    let dependent_subgraph_label = dependent_subgraph.string_label();
    let source_lmbs = active_components
        .iter()
        .map(|(_, _, lmb)| lmb)
        .collect::<Vec<_>>();
    let mut retained_loop_edges = source_lmbs
        .iter()
        .flat_map(|lmb| lmb.loop_edges.iter().copied())
        .collect::<Vec<_>>();
    retained_loop_edges.sort();
    retained_loop_edges.dedup();
    // Resolve color contractions before opening shorthand traces. Repeated
    // color spectators then share their reduced form instead of acquiring
    // distinct dummy indices during Lorentz preparation. Keep the kinematic
    // coefficients factorized throughout this independent color pass.
    let mut integrand_vakint = GS
        .erase_uv_momentum_provenance(integrand)
        .simplify_color_with(ColorSimplifySettings {
            simplify_non_color: false,
            ..Default::default()
        })
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

    // The denominator-to-propagator replacements below strip the momentum
    // wrapper without distributing the numerator. The former standalone
    // denominator rewrite was:
    // .replace(function!(
    //     GS.den,
    //     W_.prop_,
    //     function!(GS.emr_mom, W_.prop_, W_.mom_),
    //     W_.x__
    // ))
    // .with(function!(GS.den, W_.prop_, W_.mom_, W_.x__))

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

    // Separate denominator monomials before contracting their graph incidence.
    // In A*B*(c*B+d*B^2), fusing only the outside A*B would leave the
    // inner B attached to nodes that have already been contracted. Collect
    // only complete denominator atoms; numerator sums stay factored coefficients.
    let denominator_pattern = function!(GS.den, W_.x___).to_pattern();
    let denominators = integrand_vakint
        .pattern_match(&denominator_pattern, None, None)
        .map(|matched| denominator_pattern.replace_wildcards(&matched).unwrap())
        .collect::<std::collections::HashSet<_>>();
    integrand_vakint = integrand_vakint
        .coefficient_list::<i32>(&denominators.into_iter().collect::<Vec<_>>())
        .into_iter()
        .map(|(denominator, numerator)| denominator * numerator)
        .sum();
    let mut propagator_id = 1;
    let mut propagator_replacements = Vec::new();
    let mut incidence_replacements = Vec::new();
    let mut assigned_edges = std::collections::HashSet::new();

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

    for (owners, source_scope, _) in active_components {
        // Shrink only this component's omitted prefix. Active child propagators
        // retain their own incidence; its quotient attaches at the contracted
        // vertex. Contract disjoint prefixes separately, preserving their loop
        // domains instead of identifying unrelated component vertices.
        let mut contracted_nodes = std::collections::HashMap::new();
        for prefix in graph.connected_components(&source_scope.subtract(owners)) {
            let mut nodes = graph.iter_nodes_of(&prefix).map(|(node, _, _)| node);
            if let Some(first) = nodes.next() {
                for node in nodes {
                    contracted_nodes.insert(node, first);
                }
            }
        }
        for (pair, index, _data) in graph.iter_edges_of(owners) {
            if let HedgePair::Paired { source, sink } = pair {
                if !assigned_edges.insert(index) {
                    return Err(eyre!(
                        "Vakint denominator owner {index} belongs to multiple active components"
                    ));
                }
                // The former global contraction selected `first` from all
                // dependent nodes. Resolve each endpoint with this owner's
                // prefix instead, leaving active child endpoints unchanged.
                let source = graph.node_id(source);
                let sink = graph.node_id(sink);
                let contracted_source = contracted_nodes.get(&source).copied().unwrap_or(source);
                let contracted_sink = contracted_nodes.get(&sink).copied().unwrap_or(sink);
                propagator_replacements.push(Replacement::new(
                    function!(
                        GS.den,
                        usize::from(index) as i64,
                        W_.mom_,
                        W_.mass_,
                        W_.x___
                    ),
                    function!(
                        vk_prop,
                        propagator_id,
                        function!(vk_edge, usize::from(source), usize::from(sink)),
                        W_.mom_,
                        if substitute_masses_to_m_uv {
                            Atom::var(GS.m_uv_vacuum).pow(2)
                        } else {
                            Atom::var(W_.mass_)
                        },
                        1
                    ),
                ));
                incidence_replacements.push(Replacement::new(
                    function!(
                        vk_prop,
                        propagator_id,
                        function!(vk_edge, usize::from(source), usize::from(sink)),
                        W_.x___
                    ),
                    function!(
                        vk_prop,
                        propagator_id,
                        function!(
                            vk_edge,
                            usize::from(contracted_source),
                            usize::from(contracted_sink)
                        ),
                        W_.x___
                    ),
                ));
                propagator_id += 1;
            }
        }
    }
    // Edge IDs make these denominator replacements disjoint. Convert all
    // propagators together, then encode their powers once for the whole atom.
    integrand_vakint = integrand_vakint
        .replace_multiple(&propagator_replacements)
        .replace(function!(vk_prop, W_.x___, 1).pow(Atom::var(W_.e_)))
        .with(function!(vk_prop, W_.x___, -Atom::var(W_.e_)));
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_den_to_prop",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = integrand_vakint,
        "Vakint trace after denominator-to-propagator conversion"
    );

    debug_tags!(#uv, #integrated, #vakint, #graph, #dump;
        reduced = %graph.dot(dependent_subgraph),
        "Shrinking each component's omitted prefix for vakint"
    );
    // Shrink vertices only on propagators belonging to that quotient.
    integrand_vakint = integrand_vakint.replace_multiple(&incidence_replacements);
    if integrand_vakint.contains_symbol(GS.den) {
        return Err(eyre!(
            "Vakint component {reduced_label} contains a denominator without an active owner"
        ));
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
        //prop(<id>, edge(<node_source>,<node_sink>), <mom>, <mass_squared>, <power>)
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
            mass_squared: Atom,
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
                    mass_squared: m[&W_.d_].clone(),
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
            let mut mass_squared = None;
            for (_, _, d) in graph.iter_edges_of(c) {
                count += 1;

                if let Some(m) = &mass_squared
                    && m != &d.data.mass_squared
                {
                    return false;
                } else {
                    mass_squared = Some(d.data.mass_squared.clone());
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
                    .mass_squared
                    .replace(GS.m_uv_expansion)
                    .with(GS.m_uv_vacuum),
                e.data.power
            )
        }

        // println!("{}->{}", t.integral, new_integral);
        t.integral = function!(vakint::symbols::S.topo, new_integral);
        let nloops = graph.cyclotomatic_number(&graph.full_filter());
        // A cancelled denominator can remove an integration direction without
        // leaving that coordinate in the numerator. Never silently replace the
        // retained active measure by a lower-loop vacuum integral. Redundant
        // EMR variables remain allowed when they span the full retained domain.
        if nloops != retained_loop_edges.len() {
            return Err(eyre!(
                "Vakint term {term_index} for component {reduced_label} after prefix {dependent_subgraph_label} has {nloops} loop generators, expected {} retained active generators {:?}; a changed integration domain requires an explicit scaleless or factorization certificate",
                retained_loop_edges.len(),
                retained_loop_edges,
            ));
        }
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
                        // An affine vacuum routing can contain fixed source
                        // external momenta. Solve only for hard coordinates;
                        // otherwise H=Q-P can solve for P and leave Q in the
                        // integrated numerator as a spurious external variable.
                        // A child crown carrier can become an active loop of
                        // its parent. Retained hard coordinates take precedence
                        // over an external classification in another frame.
                        if usize::try_from(m[&W_.a_].as_view()).is_ok_and(|edge| {
                            let edge = EdgeIndex(edge);
                            source_lmbs.iter().any(|lmb| lmb.ext_from(edge).is_some())
                                && !source_lmbs.iter().any(|lmb| lmb.loop_edges.contains(&edge))
                        }) {
                            return;
                        }
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
        momentum_solution.ensure_free_variables_eliminated(
            term_index,
            "numerator",
            &t.numerator,
            &add_additional_args,
        )?;
        let retained_hard_coordinates = retained_loop_edges
            .iter()
            .map(|edge| function!(GS.emr_mom, usize::from(*edge)))
            .collect::<Vec<_>>();
        momentum_solution.ensure_variables_eliminated(
            term_index,
            "numerator",
            &t.numerator,
            &retained_hard_coordinates,
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
        t.numerator = Integrated::to_vakint_numerator(&t.numerator);
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

        let solutions = Atom::solve(system)
            .wrt_with_exponent::<u8, _>(variables)
            .map_err(|source| eyre!("{source}"))?;
        let [solution] = solutions.iter().as_slice() else {
            return Err(eyre!(
                "expected one Vakint momentum solution, got {} branches",
                solutions.len()
            ));
        };
        let solution = variables
            .iter()
            .map(|variable| {
                let polynomial_variable =
                    PolyVariable::try_from(variable.clone()).map_err(|source| eyre!("{source}"))?;
                if solution.free_variables().contains(&polynomial_variable) {
                    Ok(variable.clone())
                } else {
                    solution.get(&polynomial_variable).cloned().ok_or_else(|| {
                        eyre!("Vakint momentum solution has no value for {polynomial_variable}")
                    })
                }
            })
            .collect::<Result<Vec<_>>>()?;
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
    use crate::{initialisation::test_initialise, utils::symbolica_ext::Q_I};

    use super::*;

    #[test]
    fn analytic_spin_simplification_preserves_shared_color_factor() {
        use idenso::{bis, coad, cof, color_f, color_idx, color_t, gamma};
        use spenso::{mink, trace, trace_sym};

        test_initialise().unwrap();
        let a = coad!(8, color_a);
        let b = coad!(8, color_b);
        let c = coad!(8, color_c);
        let color = trace_sym!(cof!(3), color_t!(&a), color_t!(&b), color_t!(&c))
            + Atom::num(2) * color_f!(&a, &b, &c) * color_idx!(2, cof!(3));
        let spin = trace!(
            bis!(4),
            gamma!(mink!(4, mu)),
            gamma!(mink!(4, nu)),
            gamma!(mink!(4, rho)),
            gamma!(mink!(4, sigma))
        );
        let reduced_spin = simplify(&spin).unwrap();
        assert!(!reduced_spin.is_zero());
        assert!(!reduced_spin.contains_symbol(AGS.gamma));
        // Compare the literal factorization after the complete analytic spin
        // pipeline, without expanding the numerator for the assertion.
        assert_eq!(simplify(&(&color * spin)).unwrap(), color * reduced_spin);
    }

    #[test]
    fn analytic_spin_expansion_preserves_powered_scalar_products() {
        test_initialise().unwrap();
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let bis = Bispinor {}.new_rep(4);
        let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let p = function!(GS.loop_mom, 1, mink.to_symbolic([]));
        let k = function!(GS.loop_mom, 2, mink.to_symbolic([]));
        let spin = function!(AGS.gamma, &a, &b, &mu) * function!(AGS.gamma, &b, &a, &nu);
        for product in [spenso::dot!(&p, &k), spenso::g!(&p, &k)] {
            for power in [2, 3] {
                // A closed Dirac loop makes the analytic spin-expansion path
                // active. Its scalar spectator must retain its own contractions.
                let expected =
                    Atom::num(4) * spenso::g!(&mu, &nu) * spenso::dot!(&p, &k).pow(power);
                assert_eq!(
                    simplify(&(product.pow(power) * &spin))
                        .unwrap()
                        .normalize_dots(),
                    expected
                );
                // These compact metrics are tensor contractions, so they must
                // remain visible even beside the protected scalar spectator.
                assert_eq!(
                    simplify(&(product.pow(power) * &spin * spenso::g!(&mu, &p)))
                        .unwrap()
                        .normalize_dots(),
                    Atom::num(4) * function!(GS.loop_mom, 1, &nu) * spenso::dot!(&p, &k).pow(power)
                );
                assert_eq!(
                    simplify(
                        &(product.pow(power) * &spin * spenso::g!(&mu, &p) * spenso::g!(&nu, &k)),
                    )
                    .unwrap()
                    .normalize_dots(),
                    Atom::num(4) * spenso::dot!(&p, &k).pow(power + 1)
                );
            }
        }
    }

    #[test]
    fn projected_dirac_algebra_precedes_single_and_double_pole_expansion() {
        test_initialise().unwrap();
        let vakint = Vakint::new().unwrap();
        let epsilon = Atom::var(GS.dim_epsilon);
        let radial = vakint::symbols::S.dot(
            function!(vakint::symbols::S.k, 1),
            function!(vakint::symbols::S.k, 2),
        );
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let bis = Bispinor {}.new_rep(4);
        let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let c = bis.to_symbolic([Aind::new_dummy().to_atom()]);

        for project_onto_tensor_integrals in [true, false] {
            for use_dot_product_notation in [false, true] {
                let settings = vakint::VakintSettings {
                    use_dot_product_notation,
                    project_onto_tensor_integrals,
                    ..VakintSettings::default().true_settings()
                };
                let integrated = Integrated::new(&vakint, &settings);
                for closed in [false, true] {
                    // Distinct loop slashes have no repeated gamma index before
                    // angular projection. Its projector creates the contraction.
                    let numerator = function!(GS.loop_mom, 1, &mu)
                        * function!(GS.loop_mom, 2, &nu)
                        * function!(AGS.gamma, &a, &b, &mu)
                        * function!(AGS.gamma, &b, if closed { &a } else { &c }, &nu);
                    let mut term = vakint::VakintTerm {
                        integral: vakint::vakint_parse!("topo(I2L(mUVsq,1,1,1))").unwrap(),
                        numerator: Integrated::to_vakint_numerator(&numerator),
                        vectors: vec![("k".into(), 1), ("k".into(), 2)],
                    };
                    term.tensor_reduce(&vakint, &settings).unwrap();
                    let spin = if closed {
                        Atom::num(4)
                    } else {
                        function!(ETS.metric, &a, &c)
                    };
                    let late_four_dimensional = simplify(
                        &Integrated::restore_numerator(
                            Vakint::convert_to_dot_notation(&settings, term.numerator.as_view())
                                .unwrap(),
                            0,
                        )
                        .replace(GS.dim)
                        .with(4),
                    )
                    .unwrap();
                    assert!(
                        !truncate(
                            &series(
                                &(late_four_dimensional
                                    .replace(radial.to_pattern())
                                    .with(epsilon.pow(-1))
                                    - &spin * epsilon.pow(-1)),
                                1,
                            )
                            .unwrap(),
                            true,
                        )
                        .is_zero(),
                        "the old late-4D contraction must miss a nonzero finite term"
                    );
                    for evanescent_power in 0..=2 {
                        let projected = &term.numerator
                            * (Atom::var(GS.dim) - Atom::num(4)).pow(evanescent_power);
                        let resolved = integrated
                            .simplify_projected_numerator(&projected, 0)
                            .unwrap();
                        assert!(!resolved.contains_symbol(AGS.gamma));
                        assert!(!resolved.contains_symbol(SPENSO_TAG.trace));
                        assert!(!resolved.contains_symbol(GS.dim));
                        let coefficient = &spin * (Atom::num(-2) * &epsilon).pow(evanescent_power);
                        for pole in [
                            epsilon.pow(-1) + Atom::num(5) + Atom::num(7) * &epsilon,
                            Atom::num(2) * epsilon.pow(-2)
                                + Atom::num(3) * epsilon.pow(-1)
                                + Atom::num(5)
                                + Atom::num(7) * &epsilon,
                        ] {
                            // Mock only the scalar master value: the angular
                            // projection and d-dimensional Dirac algebra above use
                            // the production pipeline. This isolates the pole ×
                            // evanescent finite terms from master normalization.
                            let evaluated = resolved
                                .replace(radial.to_pattern())
                                .with(pole.to_pattern());
                            assert!(
                                series(&(evaluated - &coefficient * pole), 2)
                                    .unwrap()
                                    .to_atom()
                                    .is_zero(),
                                "project={project_onto_tensor_integrals}, closed={closed}, dot={use_dot_product_notation}, evanescent power={evanescent_power}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn projected_dirac_numerator_matches_scalar_master_before_backend_truncation() {
        test_initialise().unwrap();
        let vakint = Vakint::new().unwrap();
        let mut results = Vec::new();
        for project_onto_tensor_integrals in [true, false] {
            let settings = vakint::VakintSettings {
                // A two-loop vacuum master is retained only through its finite
                // term, so omitting the numerator's evanescent factors fails here.
                number_of_terms_in_epsilon_expansion: 3,
                project_onto_tensor_integrals,
                ..VakintSettings::default().true_settings()
            };
            let integrated = Integrated::new(&vakint, &settings);
            let mink = Minkowski {}.new_rep(GS.dim);
            let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
            let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
            let bis = Bispinor {}.new_rep(4);
            let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
            let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
            let numerator = function!(GS.loop_mom, 1, &mu)
                * function!(GS.loop_mom, 2, &nu)
                * function!(AGS.gamma, &a, &b, &mu)
                * function!(AGS.gamma, &b, &a, &nu);
            let mut projected = vakint::VakintTerm {
                integral: vakint::vakint_parse!("topo(I2L(mUVsq,1,1,1))").unwrap(),
                numerator: Integrated::to_vakint_numerator(&numerator),
                vectors: vec![("k".into(), 1), ("k".into(), 2)],
            };
            projected.tensor_reduce(&vakint, &settings).unwrap();
            let d_minus_four = Atom::var(GS.dim) - Atom::num(4);
            projected.numerator = integrated
                .simplify_projected_numerator(
                    &(&projected.numerator * (&d_minus_four + d_minus_four.pow(2))),
                    0,
                )
                .unwrap();
            let epsilon = Atom::var(GS.dim_epsilon);
            let mut scalar = vakint::VakintTerm {
                numerator: Atom::num(4)
                    * (Atom::num(-2) * &epsilon + Atom::num(4) * epsilon.pow(2))
                    * vakint::symbols::S.dot(
                        function!(vakint::symbols::S.k, 1),
                        function!(vakint::symbols::S.k, 2),
                    ),
                ..projected.clone()
            };
            let resolved_input = projected.numerator.clone();
            // The projector contains 1/(4-2ε) before symbolic-d gamma
            // contraction. Cancel the analytic oracle's rational denominators
            // exactly, using the unfactorized polynomial representation.
            assert!(
                (&resolved_input - &scalar.numerator)
                    .try_to_rational_polynomial::<_, _, u32>(&*Q_I, &*Q_I, None)
                    .expect("the analytic scalar input difference is rational over Q(i)")
                    .numerator
                    .is_zero(),
                "project={project_onto_tensor_integrals}: resolved input {resolved_input} differs from scalar oracle {}",
                scalar.numerator
            );
            projected.evaluate_integral(&vakint, &settings).unwrap();
            scalar.evaluate_integral(&vakint, &settings).unwrap();
            assert!(
                !truncate(&series(&scalar.numerator, 1).unwrap(), true).is_zero(),
                "the independently integrated evanescent scalar master has a nonzero finite term"
            );
            let difference = series(&(&projected.numerator - &scalar.numerator), 1)
                .unwrap()
                .to_atom();
            // These are analytically integrated subgraph coefficients. Expand
            // their difference to compare exact Laurent coefficients, including
            // algebraically equal backend output with different factorization.
            assert!(
                difference.expand().is_zero(),
                "project={project_onto_tensor_integrals}: resolved input {resolved_input}; projected result {}; scalar result {}; difference {difference}; expanded difference {}",
                projected.numerator,
                scalar.numerator,
                difference.expand()
            );
            results.push(projected.numerator);
        }
        assert!(
            series(&(&results[0] - &results[1]), 1)
                .unwrap()
                .to_atom()
                .expand()
                .is_zero()
        );
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
    }

    #[test]
    fn analytic_uv_spin_algebra_rejects_unsupported_residual_contractions() {
        test_initialise().unwrap();
        let bis = Bispinor {}.new_rep(4);
        let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let mu = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::new_dummy().to_atom()]);
        let unknown = function!(
            symbolica::symbol!("unsupported_dirac_factor"),
            SPENSO_TAG.chain_in,
            SPENSO_TAG.chain_out
        );
        let gamma = function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, &mu);
        for (input, message) in [
            (function!(AGS.sigma, &a, &b, &mu, &nu), "Sigma tensors"),
            (
                function!(SPENSO_TAG.trace, bis.to_symbolic([]), &unknown),
                "unresolved Dirac trace",
            ),
            (
                function!(
                    SPENSO_TAG.trace,
                    Minkowski {}.new_rep(GS.dim).to_symbolic([]),
                    &unknown
                ),
                "unresolved Lorentz trace",
            ),
        ] {
            assert!(simplify(&input).unwrap_err().to_string().contains(message));
        }
        let residual = function!(SPENSO_TAG.chain, &a, &b, &gamma, &unknown, &gamma);
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(&simplify(&residual).unwrap())
                .is_err()
        );
        // An open gamma is a valid d-dimensional tensor basis element.
        let open = simplify(&function!(AGS.gamma, &a, &b, &mu)).unwrap();
        assert!(Integrated::ensure_resolved_lorentz_contractions(&open).is_ok());
    }

    #[test]
    fn residual_tensor_contractions_are_checked_per_analytic_branch() {
        test_initialise().unwrap();
        let bis = Bispinor {}.new_rep(4);
        let [a, b, c, d] =
            std::array::from_fn::<_, 4, _>(|_| bis.to_symbolic([Aind::new_dummy().to_atom()]));
        let mink = Minkowski {}.new_rep(GS.dim);
        let [mu, nu, rho] =
            std::array::from_fn::<_, 3, _>(|_| mink.to_symbolic([Aind::new_dummy().to_atom()]));
        let left = function!(AGS.gamma, &a, &b, &mu);
        let right = function!(AGS.gamma, &c, &d, &mu);
        assert!(Integrated::ensure_resolved_lorentz_contractions(&(&left * &right)).is_err());
        let compact = function!(
            SPENSO_TAG.dot,
            function!(AGS.gamma, &a, &b, mink.to_symbolic([])),
            function!(AGS.gamma, &c, &d, mink.to_symbolic([]))
        );
        assert!(Integrated::ensure_resolved_lorentz_contractions(&compact).is_err());
        let same_chain = function!(
            SPENSO_TAG.dot,
            function!(AGS.gamma, &a, &b, mink.to_symbolic([])),
            function!(AGS.gamma, &b, &c, mink.to_symbolic([]))
        );
        let same_chain = simplify(&same_chain).unwrap();
        assert!(
            (&same_chain - Atom::var(GS.dim) * function!(ETS.metric, &a, &c))
                .expand()
                .is_zero(),
            "a compact contraction within one gamma chain must close to d times its spinor identity: {same_chain}"
        );
        assert!(Integrated::ensure_resolved_lorentz_contractions(&same_chain).is_ok());
        let free =
            &left * function!(AGS.gamma, &c, &d, &nu) + function!(AGS.gamma, &a, &b, &nu) * &right;
        assert!(Integrated::ensure_resolved_lorentz_contractions(&free).is_ok());
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(&function!(ETS.metric, &mu, &nu))
                .is_ok()
        );
        let ambiguous = function!(
            SPENSO_TAG.dot,
            function!(
                AGS.sigma,
                &a,
                &b,
                mink.to_symbolic([]),
                mink.to_symbolic([])
            ),
            function!(
                AGS.sigma,
                &c,
                &d,
                mink.to_symbolic([]),
                mink.to_symbolic([])
            )
        );
        assert!(Integrated::ensure_resolved_lorentz_contractions(&ambiguous).is_err());
        let slash = function!(
            AGS.gamma,
            &a,
            &b,
            function!(GS.emr_mom, 0, mink.to_symbolic([]))
        );
        assert!(Integrated::ensure_resolved_lorentz_contractions(&slash).is_ok());
        let slash_product = &slash
            * function!(
                AGS.gamma,
                &c,
                &d,
                function!(GS.emr_mom, 1, mink.to_symbolic([]))
            );
        assert!(Integrated::ensure_resolved_lorentz_contractions(&slash_product).is_ok());
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &(&left / (Atom::var(GS.dim) - Atom::one()))
            )
            .is_ok()
        );

        // Diagnostic tensor operator, not a physical graph/locality oracle:
        // crossed spinor closure of two contracted triple-gamma chains gives
        // 4d(-d²+6d-4). Its evanescent term changes the finite part of a pole.
        let factors = [&mu, &nu, &rho]
            .map(|index| function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, index));
        let operator = function!(
            SPENSO_TAG.chain,
            &a,
            &b,
            &factors[0],
            &factors[1],
            &factors[2]
        ) * function!(
            SPENSO_TAG.chain,
            &c,
            &d,
            &factors[0],
            &factors[1],
            &factors[2]
        );
        assert!(Integrated::ensure_resolved_lorentz_contractions(&operator).is_err());
        let closed = function!(
            SPENSO_TAG.trace,
            bis.to_symbolic([]),
            &factors[0],
            &factors[1],
            &factors[2],
            &factors[0],
            &factors[1],
            &factors[2]
        );
        let closed = Integrated::dimensionally_regularized(&simplify(&closed).unwrap());
        let epsilon = Atom::var(GS.dim_epsilon);
        assert_eq!(
            truncate(&series(&(closed / epsilon), 0).unwrap(), true),
            Atom::num(32)
        );
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
    fn finite_projection_retains_child_epsilon_terms_until_parent_poles_multiply() {
        test_initialise().unwrap();
        let epsilon = Atom::var(GS.dim_epsilon);
        let child = IntegratedCts {
            expansion: series(
                &(Atom::num(2) + Atom::num(3) * &epsilon + Atom::num(5) * epsilon.pow(2)),
                3,
            )
            .unwrap(),
            scale_power: 4,
        };
        let parent_poles = epsilon.pow(-2) + Atom::num(7) * epsilon.pow(-1);
        let nested = IntegratedCts {
            expansion: series(
                &(child.physical_finite_counterterm_atom() * parent_poles),
                1,
            )
            .unwrap(),
            scale_power: 8,
        };
        // The finite projection retains nonnegative powers for further nesting.
        let finite = nested.physical_finite_counterterm_atom();
        assert_eq!(
            finite.replace(GS.dim_epsilon).with(Atom::Zero),
            Atom::num(26)
        );
        assert!(
            (finite - Atom::num(26) - Atom::num(35) * &epsilon)
                .expand()
                .is_zero()
        );
        assert_eq!(
            nested.physical_pole_atom(),
            -Atom::num(2) * epsilon.pow(-2) - Atom::num(17) * epsilon.pow(-1),
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
    //     let minkowski = Minkowski {}.new_rep(4).to_symbolic([]);
    //     let minkowski_norm = function!(
    //         SPENSO_TAG.dot,
    //         GS.emr_vec(edge, minkowski.as_view()),
    //         GS.emr_vec(edge, minkowski.as_view())
    //     );

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
    fn vakint_preparation_reduces_color_before_opening_traces() {
        use crate::{
            dot,
            graph::{Graph, parse::IntoGraph},
        };
        use idenso::{coad, cof, color_cas, color_f, color_idx, color_t};
        use spenso::{g, trace_sym};

        test_initialise().unwrap();
        let graph: Graph = dot!(digraph color_vakint_tadpole {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0]
            a -> a [id=1 lmb_id=0]
            a -> outgoing [id=2]
        })
        .unwrap();
        let [a, b, c, d] = [
            "color_vakint_a",
            "color_vakint_b",
            "color_vakint_c",
            "color_vakint_d",
        ]
        .map(|name| coad!(8, Atom::var(symbol!(name))));
        let trace = trace_sym!(cof!(3), color_t!(&a), color_t!(&b), color_t!(&c));
        let index = color_idx!(2, cof!(3));
        let color = color_f!(&d, &b, &c) * (trace + Atom::i() / 2 * color_f!(&a, &b, &c) * &index);
        let coefficient = parse!("(color_x+color_y)^3*(color_u+color_v)");
        let denominator = function!(
            GS.den,
            1,
            function!(GS.emr_mom, 1),
            Atom::var(GS.m_uv_vacuum).pow(2)
        );
        let actual = to_vakint_integrand(
            &(color * &coefficient / denominator.pow(2)),
            &graph,
            &graph.full_filter(),
            &graph.empty_subgraph::<SuBitGraph>(),
            &[(
                graph.full_filter(),
                graph.full_filter(),
                graph.loop_momentum_basis.clone(),
            )],
            &VakintSettings {
                additional_normalization: "1".to_string(),
                ..Default::default()
            },
            true,
        )
        .unwrap();

        // The symmetric trace annihilates against f; the remaining f*f
        // contraction is C_A times the adjoint metric, independent of kinematics.
        let expected = Atom::i() / 2 * color_cas!(2, coad!(8)) * index * g!(&a, &d) * coefficient;
        assert_eq!(actual.0.len(), 1);
        assert_eq!(actual.0[0].numerator, expected);
    }

    #[test]
    fn vakint_conversion_preserves_squared_masses() {
        use crate::{
            dot,
            graph::{Graph, parse::IntoGraph},
        };

        test_initialise().unwrap();
        let graph: Graph = dot!(digraph vakint_mass_tadpole {
            edge [num=1 mass=2]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0]
            a -> a [id=1 lmb_id=0]
            a -> outgoing [id=2]
        })
        .unwrap();
        let momentum = function!(GS.emr_mom, 1);
        let indexed_momentum =
            function!(GS.emr_mom, 1, Minkowski {}.new_rep(GS.dim).to_symbolic([]));
        let settings = VakintSettings {
            additional_normalization: "1".to_string(),
            ..Default::default()
        };
        // Both wrappers carry squared masses. A nonunit physical mass detects
        // an accidental second square, while the UV substitution remains MUV^2.
        for mass in [Atom::num(2), parse!("m_phys")] {
            let mass_squared = mass.pow(2);
            let denominator = function!(
                GS.den,
                1,
                &momentum,
                &mass_squared,
                function!(SPENSO_TAG.dot, &indexed_momentum, &indexed_momentum) - &mass_squared
            );
            for substitute_masses_to_m_uv in [false, true] {
                let actual = to_vakint_integrand(
                    &denominator.pow(-3),
                    &graph,
                    &graph.full_filter(),
                    &graph.empty_subgraph::<SuBitGraph>(),
                    &[(
                        graph.full_filter(),
                        graph.full_filter(),
                        graph.loop_momentum_basis.clone(),
                    )],
                    &settings,
                    substitute_masses_to_m_uv,
                )
                .unwrap();
                assert_eq!(actual.0.len(), 1);
                assert_eq!(actual.0[0].numerator, Atom::one());
                assert_eq!(
                    actual.0[0].integral,
                    function!(
                        vakint::symbols::S.topo,
                        function!(
                            vakint::symbols::S.prop,
                            1,
                            function!(vakint::symbols::S.edge, 0, 0),
                            function!(vakint::symbols::S.k, 0),
                            if substitute_masses_to_m_uv {
                                Atom::var(GS.m_uv_vacuum).pow(2)
                            } else {
                                mass_squared.clone()
                            },
                            3
                        )
                    )
                );
            }
        }
    }

    #[test]
    fn vakint_conversion_rejects_a_lost_active_loop() {
        use crate::{
            dot,
            graph::{Graph, parse::IntoGraph},
        };

        test_initialise().unwrap();
        let graph: Graph = dot!(digraph cancelled_tadpole_loop {
            edge [num=1 mass=2]
            node [num=1]
            a -> a [id=0 lmb_id=0]
            a -> a [id=1 lmb_id=1]
        })
        .unwrap();
        let remaining = graph.get_edge_subgraph(EdgeIndex(1));
        // This is the denominator domain after the other tadpole propagator
        // cancels against its edge-local numerator. The original two-loop
        // measure cannot be replaced by the remaining massive one-loop value.
        let input = graph.denominator(&remaining, |_| -3);
        let error = to_vakint_integrand(
            &input,
            &graph,
            &graph.full_filter(),
            &graph.empty_subgraph::<SuBitGraph>(),
            &[(
                graph.full_filter(),
                graph.full_filter(),
                graph.loop_momentum_basis.clone(),
            )],
            &VakintSettings::default(),
            false,
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("has 1 loop generators, expected 2")
        );
        let remaining_lmb = graph.lmb_of(&remaining);
        let converted = to_vakint_integrand(
            &input,
            &graph,
            &graph.full_filter(),
            &graph.empty_subgraph::<SuBitGraph>(),
            &[(
                graph.full_filter(),
                graph.full_filter(),
                remaining_lmb.clone(),
            )],
            &VakintSettings::default(),
            false,
        )
        .unwrap();
        assert_eq!(converted.0.len(), 1);
    }

    #[test]
    fn affine_vakint_routing_keeps_external_momenta_fixed() {
        use crate::{
            dot,
            graph::{Graph, parse::IntoGraph},
        };

        test_initialise().unwrap();
        let graph: Graph = dot!(digraph affine_vakint_tadpole {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0]
            a -> a [id=1 lmb_id=0]
            a -> outgoing [id=2]
        })
        .unwrap();
        let hard = function!(GS.emr_mom, 1) - function!(GS.emr_mom, 0);
        let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let indexed_hard =
            function!(GS.emr_mom, 1, &minkowski) - function!(GS.emr_mom, 0, &minkowski);
        let numerator = function!(
            SPENSO_TAG.dot,
            function!(GS.emr_mom, 0, &minkowski),
            &indexed_hard
        ) * function!(
            SPENSO_TAG.dot,
            function!(GS.emr_mom, 2, &minkowski),
            &indexed_hard
        );
        let settings = VakintSettings {
            additional_normalization: "1".to_string(),
            ..Default::default()
        };
        // A retained parent loop overrides an earlier component's external
        // classification. This adapter diagnostic changes coordinate metadata
        // only; both signs must still preserve the fixed physical externals.
        let mut child_frame = graph.empty_lmb();
        child_frame.ext_edges.push(EdgeIndex(1));
        for source_lmbs in [
            vec![&graph.loop_momentum_basis],
            vec![&child_frame, &graph.loop_momentum_basis],
        ] {
            for sign in [-1, 1] {
                let denominator = function!(
                    GS.den,
                    1,
                    Atom::num(sign) * &hard,
                    Atom::var(GS.m_uv_vacuum).pow(2),
                    function!(SPENSO_TAG.dot, &indexed_hard, &indexed_hard)
                        - Atom::var(GS.m_uv_vacuum).pow(2)
                );
                let actual = to_vakint_integrand(
                    &(&numerator / denominator.pow(3)),
                    &graph,
                    &graph.full_filter(),
                    &graph.empty_subgraph::<SuBitGraph>(),
                    &source_lmbs
                        .iter()
                        .map(|lmb| {
                            let owners = if lmb.loop_edges.is_empty() {
                                graph.empty_subgraph::<SuBitGraph>()
                            } else {
                                graph.full_filter()
                            };
                            (owners.clone(), owners, (*lmb).clone())
                        })
                        .collect::<Vec<_>>(),
                    &settings,
                    true,
                )
                .unwrap();
                // Compare the complete rational integrand through Vakint's public
                // expression conversion, independent of topology IDs and term layout.
                let actual = Atom::from(actual)
                    .replace(function!(vakint::symbols::S.topo, W_.x_))
                    .with(W_.x_)
                    .replace(function!(
                        vakint::symbols::S.prop,
                        W_.a_,
                        W_.b_,
                        W_.c_,
                        W_.d_,
                        W_.e_
                    ))
                    .with(
                        (function!(vakint::symbols::S.dot, W_.c_, W_.c_) - Atom::var(W_.d_))
                            .pow(-Atom::var(W_.e_)),
                    );
                let expected_numerator = function!(
                    vakint::symbols::S.dot,
                    function!(vakint::symbols::S.p, 0),
                    function!(vakint::symbols::S.k, 0)
                ) * function!(
                    vakint::symbols::S.dot,
                    function!(vakint::symbols::S.p, 2),
                    function!(vakint::symbols::S.k, 0)
                );
                let radial = function!(
                    vakint::symbols::S.dot,
                    function!(vakint::symbols::S.k, 0),
                    function!(vakint::symbols::S.k, 0)
                ) - Atom::var(GS.m_uv_vacuum).pow(2);
                assert!(
                    (actual.collect_factors()
                        - (expected_numerator / radial.pow(3)).collect_factors())
                    .is_zero(),
                    "the complete integrand must retain the fixed external momenta for either D(H) spelling"
                );
            }
        }
    }

    #[test]
    fn nested_active_sunset_keeps_its_two_loop_incidence() {
        use crate::{
            dot,
            graph::{Graph, parse::IntoGraph},
            momentum::sample::LoopIndex,
        };

        test_initialise().unwrap();
        let graph: Graph = dot!(digraph nested_active_sunset {
            edge [num=1 mass=1]
            node [num=1]
            ext [style=invis]
            ext -> a:0 [id=0]
            d:1 -> ext [id=1]
            a -> b [id=2]
            b -> c [id=3 lmb_id=0]
            b -> c [id=4]
            c -> d [id=5]
            a -> d [id=6 lmb_id=2]
            b -> c [id=7 lmb_id=1]
        })
        .unwrap();
        let full = graph
            .full_filter()
            .subtract(&graph.external_filter::<SuBitGraph>());
        let mut child = graph.empty_subgraph::<SuBitGraph>();
        for edge in [3, 4, 7] {
            child.add(graph[&EdgeIndex(edge)].1);
        }
        let quotient = full.subtract(&child);
        let child_lmb = graph
            .try_compatible_sub_lmb(
                &child,
                graph.dummy_less_full_crown(&child),
                &graph.loop_momentum_basis,
            )
            .unwrap();
        let mut quotient_lmb = graph.loop_momentum_basis.clone();
        for index in (0..quotient_lmb.loop_edges.len()).rev() {
            let index = LoopIndex(index);
            if child_lmb
                .loop_edges
                .contains(&quotient_lmb.loop_edges[index])
            {
                quotient_lmb.put_loop_to_ext(index);
            }
        }
        let [k, l, p] = [3, 7, 6].map(|edge| function!(GS.emr_mom, edge));
        let mass_squared = Atom::var(GS.m_uv_vacuum).pow(2);
        let child_denominators = function!(GS.den, 3, &k, &mass_squared)
            * function!(GS.den, 7, &l, &mass_squared)
            * function!(GS.den, 4, -&k - &l, &mass_squared).pow(3);
        let quotient_denominators = function!(GS.den, 2, -&p, &mass_squared)
            * function!(GS.den, 5, -&p, &mass_squared)
            * function!(GS.den, 6, &p, &mass_squared);
        let numerator = parse!("(a+b)*(c+d)");
        for child_active in [true, false] {
            let mut components = Vec::new();
            let input = if child_active {
                components.push((child.clone(), child.clone(), child_lmb.clone()));
                &numerator / (&child_denominators * &quotient_denominators)
            } else {
                &numerator / &quotient_denominators
            };
            components.push((quotient.clone(), full.clone(), quotient_lmb.clone()));
            let converted = to_vakint_integrand(
                &input,
                &graph,
                &full,
                &child,
                &components,
                &VakintSettings::default(),
                false,
            )
            .unwrap();
            assert_eq!(converted.0.len(), 1);
            let term = &converted.0[0];
            assert_eq!(term.numerator, numerator);
            let propagator = function!(
                vakint::symbols::S.prop,
                W_.a_,
                W_.b_,
                W_.mom_,
                W_.mass_,
                W_.e_
            )
            .to_pattern();
            let mut powers = term
                .integral
                .pattern_match(&propagator, None, None)
                .map(|matched| {
                    assert_eq!(matched[&W_.mass_], mass_squared);
                    assert!(!matched[&W_.mom_].contains_symbol(vakint::symbols::S.p));
                    i32::try_from(matched[&W_.e_].as_view()).unwrap()
                })
                .collect::<Vec<_>>();
            powers.sort();
            assert_eq!(
                powers,
                if child_active {
                    vec![1, 1, 3, 3]
                } else {
                    vec![3]
                }
            );
        }
    }

    #[test]
    fn nested_vacuum_denominators_preserve_incidence_and_factorized_numerators() {
        use crate::{
            dot,
            graph::{Graph, parse::IntoGraph},
        };

        test_initialise().unwrap();
        let graph: Graph = dot!(digraph nested_vacuum_bubble {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0]
            b -> outgoing [id=1]
            a -> b [id=2 lmb_id=0]
            a -> b [id=3]
        })
        .unwrap();
        let denominators = [(2, 1), (3, -1)].map(|(edge, sign)| {
            function!(
                GS.den,
                edge,
                Atom::num(sign) * function!(GS.emr_mom, 2),
                Atom::var(GS.m_uv_vacuum).pow(2)
            )
        });
        let numerator = parse!("(a+b)*(c+d)");
        let [first, second] = [parse!("s"), parse!("t")];
        // Include powers beyond signed eight-bit range, matching the i32
        // propagator powers used by the incidence conversion below.
        for base_power in [1, 127] {
            let input = &numerator / (&denominators[0] * denominators[1].pow(base_power))
                * (&first / &denominators[1] + &second / denominators[1].pow(2));
            let actual = to_vakint_integrand(
                &input,
                &graph,
                &graph.full_filter(),
                &graph.empty_subgraph::<SuBitGraph>(),
                &[(
                    graph.full_filter(),
                    graph.full_filter(),
                    graph.loop_momentum_basis.clone(),
                )],
                &VakintSettings {
                    // This conversion oracle uses unit loop normalization in every prefix.
                    additional_normalization: "1".to_string(),
                    ..VakintSettings::default()
                },
                true,
            )
            .unwrap();
            let actual = Atom::from(actual)
                .replace(function!(vakint::symbols::S.topo, W_.x_))
                .with(W_.x_)
                .replace(function!(
                    vakint::symbols::S.prop,
                    W_.a_,
                    W_.b_,
                    W_.c_,
                    W_.d_,
                    W_.e_
                ))
                .with(
                    (function!(vakint::symbols::S.dot, W_.c_, W_.c_) - Atom::var(W_.d_))
                        .pow(-Atom::var(W_.e_)),
                );
            let radial = function!(
                vakint::symbols::S.dot,
                function!(vakint::symbols::S.k, 0),
                function!(vakint::symbols::S.k, 0)
            ) - Atom::var(GS.m_uv_vacuum).pow(2);
            let expected = &numerator
                * (&first / radial.pow(base_power + 2) + &second / radial.pow(base_power + 3));
            assert_eq!(actual.collect_factors(), expected.collect_factors());
        }
    }
}
