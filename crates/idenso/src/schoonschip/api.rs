use spenso::{
    network::{
        ExecutionResult, Sequential, TensorOrScalarOrKey,
        library::{DummyLibrary, function_lib::Wrap, symbolic::ETS},
        parsing::ParseSettings,
        tags::SPENSO_TAG as T,
    },
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::slot::{AbsInd, DummyAind, ParseableAind},
};

use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::MatchStack,
};

use crate::{
    W_,
    metric::not_slot,
    tensor::{SymbolicNetParse, SymbolicTensor},
};

use super::{
    chain_like::simplify_chain_like_metric_products,
    contraction::{
        ORDER_MIN_LARGEST_OPERAND_BYTES, ORDER_MIN_PRODUCT_BYTES, ORDER_MIN_PRODUCT_TERMS,
        ORDER_SMALLEST_DEGREE_MIN_LARGEST_OPERAND_BYTES, ORDER_SMALLEST_DEGREE_MIN_PRODUCT_BYTES,
        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_TERMS, SchoonschipExpressionOrder,
        SchoonschipLargestDegree, SchoonschipSmallestDegree,
    },
    settings::{
        SchoonschipContractionOrder, SchoonschipMode, SchoonschipSettings, SchoonschipTraversal,
    },
    utils::{TRACE_SCHOONSCHIP, trace_schoonschip_pattern_misses, trace_schoonschip_patterns},
};

fn positive_even_power(exp: AtomView<'_>) -> bool {
    matches!(i64::try_from(exp), Ok(exp) if exp > 0 && exp % 2 == 0)
}

fn positive_odd_power(exp: AtomView<'_>) -> bool {
    matches!(i64::try_from(exp), Ok(exp) if exp > 0 && exp % 2 == 1)
}

fn matched_power(matches: &MatchStack<'_>, power: symbolica::atom::Symbol) -> i64 {
    i64::try_from(&matches.get(power).unwrap().to_atom()).unwrap()
}

fn matched_pattern(pattern: &Atom, matches: &MatchStack<'_>) -> Atom {
    pattern.to_pattern().replace_wildcards_with_matches(matches)
}

fn pow_if_needed(base: Atom, exponent: i64) -> Atom {
    match exponent {
        0 => Atom::num(1),
        1 => base,
        _ => base.pow(Atom::num(exponent)),
    }
}

fn even_power_replacement(base: Atom, exponent: i64) -> Atom {
    pow_if_needed(base, exponent / 2)
}

fn odd_power_replacement(base: Atom, square: Atom, exponent: i64) -> Atom {
    let half_power = exponent / 2;
    if half_power == 0 {
        base
    } else {
        pow_if_needed(square, half_power) * base
    }
}
pub trait Schoonschip {
    fn schoonschip(&self) -> Atom;

    fn schoonschip_with_settings(&self, settings: &SchoonschipSettings) -> Atom;

    fn normalize_dots(&self) -> Atom;
    fn schoonschip_once_with_net<
        const EXPANDSUMS: bool,
        const RECURSE: bool,
        const DEPTH_FIRST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom;

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        const DEEPEST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom;
}

impl Schoonschip for Atom {
    fn schoonschip(&self) -> Atom {
        self.as_view().schoonschip()
    }

    fn schoonschip_with_settings(&self, settings: &SchoonschipSettings) -> Atom {
        self.as_view().schoonschip_with_settings(settings)
    }

    fn normalize_dots(&self) -> Atom {
        self.as_view().normalize_dots()
    }

    fn schoonschip_once_with_net<
        const EXPANDSUMS: bool,
        const RECURSE: bool,
        const DEPTH_FIRST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        self.as_view()
            .schoonschip_once_with_net::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(settings)
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        const DEEPEST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        self.as_view()
            .schoonschip_with_net::<EXPANDSUMS, DEEPEST, Aind>(settings)
    }
}

impl Schoonschip for AtomView<'_> {
    fn normalize_dots(&self) -> Atom {
        let index_cond = T.index_fiter(W_.i_);
        let other_index_cond = T.index_fiter(W_.j_);
        let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
        let stripped = T.rep_::<0, _>([W_.d_]);
        let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
        let self_dual_j = T.self_dual_::<0, _>([W_.d_, W_.j_]);
        let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
        let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
        let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);
        let dualizable_dual_j = T.dualizable_dual_::<0, _>([W_.d_, W_.j_]);
        fn rank_one_with_slot(slot: &Atom) -> Atom {
            T.rank1_::<0, _>([Atom::var(W_.c___), slot.clone()])
        }
        let self_dual_vector = rank_one_with_slot(&self_dual);
        let self_dual_vector_stripped = rank_one_with_slot(&self_dual_stripped);
        let dualizable_vector = rank_one_with_slot(&dualizable);
        let dualizable_dual_vector = rank_one_with_slot(&dualizable_dual);
        let dualizable_vector_stripped = rank_one_with_slot(&dualizable_stripped);

        let self_dual_square = ETS.metric(&self_dual_vector_stripped, &self_dual_vector_stripped);
        let self_dual_metric = function!(ETS.metric, &self_dual, &self_dual_j);
        let dualizable_metric = function!(ETS.metric, &dualizable, &dualizable_dual_j);

        // p(..,q(...,rep(d)))-> g(p(..,rep(d)),q(..,rep(d)))
        self.replace(T.rank1_::<0, _>([
            Atom::var(W_.c___),
            T.rank1_::<1, _>([&Atom::var(W_.a___), &stripped]),
        ]))
        .with(ETS.metric(
            T.rank1_::<0, _>([&Atom::var(W_.c___), &stripped]),
            T.rank1_::<1, _>([&Atom::var(W_.a___), &stripped]),
        ))
        // g(rep(d,nu),p(..,rep(d)))->p(..,rep(d,nu))
        .replace(function!(
            ETS.metric,
            &self_dual,
            &self_dual_vector_stripped
        ))
        .when(index_cond.clone())
        .repeat()
        .with(self_dual_vector.clone())
        .replace(function!(
            ETS.metric,
            &dualizable,
            &dualizable_vector_stripped
        ))
        .when(index_cond.clone())
        .repeat()
        .with(dualizable_vector.clone())
        .replace(function!(
            ETS.metric,
            &dualizable_dual,
            &dualizable_vector_stripped
        ))
        .when(index_cond.clone())
        .repeat()
        .with(dualizable_dual_vector.clone())
        // g(rep(d,nu),p(..))->p(..,rep(d,nu))
        .replace(function!(
            ETS.metric,
            &self_dual,
            T.rank1_::<0, _>([Atom::var(W_.c___)])
        ))
        .when(index_cond.clone())
        .repeat()
        .with(self_dual_vector.clone())
        // Powers of a schoonschipped vector are normalized by parity:
        // p(mu)^(2n) -> g(p,p)^n and
        // p(mu)^(2n+1) -> g(p,p)^n p(mu).
        .replace(self_dual_vector.clone().pow(Atom::var(W_.n_)))
        .when(index_cond.clone() & W_.n_.filter_single(positive_even_power))
        .with_map({
            let self_dual_square = self_dual_square.clone();
            move |matches| {
                even_power_replacement(
                    matched_pattern(&self_dual_square, matches),
                    matched_power(matches, W_.n_),
                )
            }
        })
        .replace(self_dual_vector.clone().pow(Atom::var(W_.n_)))
        .when(index_cond.clone() & W_.n_.filter_single(positive_odd_power))
        .with_map({
            let self_dual_square = self_dual_square.clone();
            let self_dual_vector = self_dual_vector.clone();
            move |matches| {
                odd_power_replacement(
                    matched_pattern(&self_dual_vector, matches),
                    matched_pattern(&self_dual_square, matches),
                    matched_power(matches, W_.n_),
                )
            }
        })
        // Metric powers follow the same parity split:
        // g(mu,nu)^(2n) -> d^n and
        // g(mu,nu)^(2n+1) -> d^n g(mu,nu).
        .replace(self_dual_metric.clone().pow(Atom::var(W_.n_)))
        .when(
            index_cond.clone()
                & other_index_cond.clone()
                & W_.n_.filter_single(positive_even_power),
        )
        .with_map(move |matches| {
            even_power_replacement(
                matched_pattern(&Atom::var(W_.d_), matches),
                matched_power(matches, W_.n_),
            )
        })
        .replace(self_dual_metric.clone().pow(Atom::var(W_.n_)))
        .when(
            index_cond.clone() & other_index_cond.clone() & W_.n_.filter_single(positive_odd_power),
        )
        .with_map({
            let self_dual_metric = self_dual_metric.clone();
            move |matches| {
                odd_power_replacement(
                    matched_pattern(&self_dual_metric, matches),
                    matched_pattern(&Atom::var(W_.d_), matches),
                    matched_power(matches, W_.n_),
                )
            }
        })
        .replace(dualizable_metric.clone().pow(Atom::var(W_.n_)))
        .when(
            index_cond.clone()
                & other_index_cond.clone()
                & W_.n_.filter_single(positive_even_power),
        )
        .with_map(move |matches| {
            even_power_replacement(
                matched_pattern(&Atom::var(W_.d_), matches),
                matched_power(matches, W_.n_),
            )
        })
        .replace(dualizable_metric.clone().pow(Atom::var(W_.n_)))
        .when(index_cond.clone() & other_index_cond & W_.n_.filter_single(positive_odd_power))
        .with_map({
            let dualizable_metric = dualizable_metric.clone();
            move |matches| {
                odd_power_replacement(
                    matched_pattern(&dualizable_metric, matches),
                    matched_pattern(&Atom::var(W_.d_), matches),
                    matched_power(matches, W_.n_),
                )
            }
        })
        // Plain metric traces are the remaining non-power case.
        .replace(function!(ETS.metric, &self_dual, &self_dual))
        .when(&index_cond)
        .with(Atom::var(W_.d_))
        .replace(function!(ETS.metric, &dualizable, &dualizable_dual))
        .when(&index_cond)
        .with(Atom::var(W_.d_))
    }

    fn schoonschip(&self) -> Atom {
        self.schoonschip_with_settings(&SchoonschipSettings::default())
    }

    fn schoonschip_with_settings(&self, settings: &SchoonschipSettings) -> Atom {
        let index_cond = T.index_fiter(W_.i_);
        let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
        let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
        let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
        let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
        let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

        let metric_self_dual = function!(ETS.metric, W_.c_, &self_dual);
        let metric_self_dual_reversed = function!(ETS.metric, &self_dual, W_.c_);
        let function_with_self_dual = function!(W_.a_, W_.a___, &self_dual, W_.b___);
        let function_with_replacement = function!(W_.a_, W_.a___, W_.c_, W_.b___);
        let metric_dualizable = function!(ETS.metric, W_.c_, &dualizable);
        let metric_dualizable_reversed = function!(ETS.metric, &dualizable, W_.c_);
        let metric_dualizable_dual = function!(ETS.metric, W_.c_, &dualizable_dual);
        let metric_dualizable_dual_reversed = function!(ETS.metric, &dualizable_dual, W_.c_);
        let function_with_dualizable = function!(W_.a_, W_.a___, &dualizable, W_.b___);
        let function_with_dualizable_dual = function!(W_.a_, W_.a___, &dualizable_dual, W_.b___);

        // The broad bare-symbolic pass may use product
        // patterns over plain functions; the network path deliberately calls
        // `normalize_dots()` instead so these broad patterns do not pre-empt
        // network-backed contractions.
        let metric_simplified = self
            .replace(metric_self_dual * function_with_self_dual.clone())
            .repeat()
            .with(function_with_replacement.clone())
            .replace(metric_self_dual_reversed * function_with_self_dual)
            .repeat()
            .with(function_with_replacement.clone())
            .replace(metric_dualizable * function_with_dualizable_dual.clone())
            .repeat()
            .with(function_with_replacement.clone())
            .replace(metric_dualizable_reversed * function_with_dualizable_dual)
            .repeat()
            .with(function_with_replacement.clone())
            .replace(metric_dualizable_dual * function_with_dualizable.clone())
            .repeat()
            .with(function_with_replacement.clone())
            .replace(metric_dualizable_dual_reversed * function_with_dualizable)
            .repeat()
            .with(function_with_replacement);

        let broad_self_dual_product_pattern =
            function!(W_.f_, W_.a___, &self_dual) * function!(W_.g_, W_.b___, &self_dual);
        let broad_self_dual_product_replacement = ETS.metric(
            function!(W_.f_, W_.a___, &self_dual_stripped),
            function!(W_.g_, W_.b___, &self_dual_stripped),
        );
        let after_broad_self_dual_product = metric_simplified
            .replace(broad_self_dual_product_pattern.clone())
            .when(index_cond.clone() & not_slot(W_.a___) & not_slot(W_.b___))
            .with(broad_self_dual_product_replacement.clone());

        let trace_patterns = trace_schoonschip_patterns();
        let trace_pattern_misses = trace_schoonschip_pattern_misses();
        if trace_patterns
            && (after_broad_self_dual_product != metric_simplified || trace_pattern_misses)
        {
            eprintln!(
                "schoonschip pattern=broad_self_dual_product changed={} pattern={} replacement={} before={} after={}",
                after_broad_self_dual_product != metric_simplified,
                broad_self_dual_product_pattern,
                broad_self_dual_product_replacement,
                metric_simplified,
                after_broad_self_dual_product
            );
        }

        let self_dual_power_pattern = function!(W_.f_, W_.a___, &self_dual).pow(Atom::num(2));
        let self_dual_power_replacement = ETS.metric(
            function!(W_.f_, W_.a___, &self_dual_stripped),
            function!(W_.f_, W_.a___, &self_dual_stripped),
        );
        let after_self_dual_power = after_broad_self_dual_product
            .replace(self_dual_power_pattern.clone())
            .when(index_cond.clone() & not_slot(W_.a___))
            .with(self_dual_power_replacement.clone());

        if trace_patterns
            && (after_self_dual_power != after_broad_self_dual_product || trace_pattern_misses)
        {
            eprintln!(
                "schoonschip pattern=self_dual_power changed={} pattern={} replacement={} before={} after={}",
                after_self_dual_power != after_broad_self_dual_product,
                self_dual_power_pattern,
                self_dual_power_replacement,
                after_broad_self_dual_product,
                after_self_dual_power
            );
        }

        let simplified = after_self_dual_power
            .replace(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual])
                    * T.rank1_::<1, _>([&Atom::var(W_.a___), &self_dual]),
            )
            .when(&index_cond)
            .with(ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                T.rank1_::<1, _>([&Atom::var(W_.a___), &self_dual_stripped]),
            ))
            .replace(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual])
                    .pow(Atom::num(2)),
            )
            .when(&index_cond)
            .with(ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
            ))
            .replace(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable])
                    * T.rank1_::<1, _>([&Atom::var(W_.a___), &dualizable_dual]),
            )
            .when(&index_cond)
            .with(ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                T.rank1_::<1, _>([&Atom::var(W_.a___), &dualizable_stripped]),
            ))
            // Bare product contraction: vector(rep(d,i)) * T(..,rep(d,i),..)
            // becomes T(..,vector(rep(d)),..). Network contraction has a more
            // structured version of this rule and does not use this pass.
            .replace(
                function!(W_.a_, W_.a___, &self_dual, W_.b___)
                    * T.rank1_::<0, _>([Atom::var(W_.c___), self_dual]),
            )
            .when(&index_cond)
            .repeat()
            .with(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                W_.b___
            ))
            .replace(
                function!(W_.a_, W_.a___, &dualizable, W_.b___)
                    * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_dual]),
            )
            .when(&index_cond)
            .repeat()
            .with(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                W_.b___
            ))
            .replace(
                function!(W_.a_, W_.a___, &dualizable_dual, W_.b___)
                    * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable]),
            )
            .repeat()
            .with(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_stripped]),
                W_.b___
            ));

        let simplified = if settings.simplify_chain_like_functions {
            simplify_chain_like_metric_products(simplified.as_view())
        } else {
            simplified
        };

        simplified.normalize_dots()
    }

    fn schoonschip_once_with_net<
        const EXPANDSUMS: bool,
        const RECURSE: bool,
        const DEPTH_FIRST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        let normalized = self.normalize_dots();
        let mut net = normalized
            .as_view()
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                depth_limit: settings.depth_limit,
                take_first_term_from_sum: false,
                parse_inner_products: settings.parse_inner_products,
                parse_composite_scalars_as_tensors: RECURSE,
                ..Default::default()
            })
            .unwrap();
        let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();

        match settings.contraction_order {
            SchoonschipContractionOrder::SmallestDegree => net
                .execute::<
                    Sequential,
                    SchoonschipSmallestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::LargestDegree => net
                .execute::<
                    Sequential,
                    SchoonschipLargestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::MinLargestOperandBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_MIN_LARGEST_OPERAND_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::MinProductTerms => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_MIN_PRODUCT_TERMS,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::MinProductBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_MIN_PRODUCT_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::SmallestDegreeMinLargestOperandBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_SMALLEST_DEGREE_MIN_LARGEST_OPERAND_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::SmallestDegreeMinProductTerms => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_TERMS,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::SmallestDegreeMinProductBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
        };

        match net.result().unwrap() {
            ExecutionResult::One => Atom::num(1),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(a) => match a {
                TensorOrScalarOrKey::Key { .. } => panic!("unexpected library key result"),
                TensorOrScalarOrKey::Scalar(s) => s.clone(),
                TensorOrScalarOrKey::Tensor { tensor, .. } => tensor.expression.clone(),
            },
        }
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        const DEEPEST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        if let AtomView::Add(add) = self {
            return add
                .iter()
                .map(|term| term.schoonschip_with_net::<EXPANDSUMS, DEEPEST, Aind>(settings))
                .fold(Atom::Zero, |sum, term| sum + term);
        }

        let new = if settings.expand_contracted_sums {
            match (settings.mode, DEEPEST) {
                (SchoonschipMode::SinglePass, _) | (_, false) => {
                    self.schoonschip_once_with_net::<true, false, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst), true) => {
                    self.schoonschip_once_with_net::<true, true, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst), true) => {
                    self.schoonschip_once_with_net::<true, true, false, Aind>(settings)
                }
            }
        } else {
            match (settings.mode, DEEPEST) {
                (SchoonschipMode::SinglePass, _) | (_, false) => {
                    self.schoonschip_once_with_net::<EXPANDSUMS, false, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst), true) => {
                    self.schoonschip_once_with_net::<EXPANDSUMS, true, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst), true) => {
                    self.schoonschip_once_with_net::<EXPANDSUMS, true, false, Aind>(settings)
                }
            }
        };

        if TRACE_SCHOONSCHIP {
            println!(
                "New: {}",
                new.printer(SpensoPrintSettings::compact().nice_symbolica())
            );
        }

        new.normalize_dots()
    }
}
