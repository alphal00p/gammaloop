use spenso::{
    network::{
        ExecutionResult, Sequential, TensorOrScalarOrKey,
        library::{DummyLibrary, function_lib::Wrap, symbolic::ETS},
        parsing::{ParseSettings, ShorthandParsing, StructureInferenceMode},
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
    tensor::{SymbolicNetParse, SymbolicTensor},
};

use super::{
    contraction::{
        ORDER_MIN_LARGEST_OPERAND_BYTES, ORDER_MIN_PRODUCT_BYTES, ORDER_MIN_PRODUCT_TERMS,
        ORDER_SMALLEST_DEGREE_MIN_LARGEST_OPERAND_BYTES, ORDER_SMALLEST_DEGREE_MIN_PRODUCT_BYTES,
        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_TERMS, SchoonschipExpressionOrder,
        SchoonschipLargestDegree, SchoonschipSmallestDegree,
    },
    settings::{
        SchoonschipContractionOrder, SchoonschipMode, SchoonschipSettings, SchoonschipTraversal,
    },
    utils::TRACE_SCHOONSCHIP,
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

    fn schoonschip_net<Aind: AbsInd + DummyAind + ParseableAind + 'static>(&self) -> Atom;

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
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

    fn schoonschip_net<Aind: AbsInd + DummyAind + ParseableAind + 'static>(&self) -> Atom {
        self.as_view().schoonschip_net::<Aind>()
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        self.as_view()
            .schoonschip_with_net::<EXPANDSUMS, Aind>(settings)
    }
}

struct NetworkSchoonschip<'a> {
    settings: &'a SchoonschipSettings,
}

impl NetworkSchoonschip<'_> {
    fn run<const EXPANDSUMS: bool, Aind>(&self, view: AtomView<'_>) -> Atom
    where
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    {
        let mut current = view.to_owned();
        loop {
            let next = self.apply::<EXPANDSUMS, Aind>(current.as_view());
            if self.settings.mode == SchoonschipMode::SinglePass || next == current {
                return next;
            }
            current = next;
        }
    }

    fn apply<const EXPANDSUMS: bool, Aind>(&self, view: AtomView<'_>) -> Atom
    where
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    {
        if let AtomView::Add(add) = view {
            return add
                .iter()
                .map(|term| self.apply::<EXPANDSUMS, Aind>(term))
                .fold(Atom::Zero, |sum, term| sum + term);
        }

        let new = match (self.settings.expand_contracted_sums, self.settings.mode) {
            (true, SchoonschipMode::SinglePass) => self.run_once::<true, false, true, Aind>(view),
            (true, SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst)) => {
                self.run_once::<true, true, true, Aind>(view)
            }
            (true, SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst)) => {
                self.run_once::<true, true, false, Aind>(view)
            }
            (false, SchoonschipMode::SinglePass) => {
                self.run_once::<EXPANDSUMS, false, true, Aind>(view)
            }
            (false, SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst)) => {
                self.run_once::<EXPANDSUMS, true, true, Aind>(view)
            }
            (false, SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst)) => {
                self.run_once::<EXPANDSUMS, true, false, Aind>(view)
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

    fn run_once<const EXPANDSUMS: bool, const RECURSE: bool, const DEPTH_FIRST: bool, Aind>(
        &self,
        view: AtomView<'_>,
    ) -> Atom
    where
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    {
        let normalized = view.normalize_dots();
        let mut net = normalized
            .as_view()
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                depth_limit: self.settings.depth_limit,
                take_first_term_from_sum: false,
                shorthand_parsing: if self.settings.parse_inner_products {
                    ShorthandParsing::Expand
                } else {
                    ShorthandParsing::Opaque {
                        inference: StructureInferenceMode::Fast,
                    }
                },
                parse_composite_scalars_as_tensors: RECURSE,
                ..Default::default()
            })
            .unwrap();
        let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();

        match self.settings.contraction_order {
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
        super::with_settings::schoonschip_with_settings(*self, settings)
    }

    fn schoonschip_net<Aind: AbsInd + DummyAind + ParseableAind + 'static>(&self) -> Atom {
        self.schoonschip_with_net::<false, Aind>(&SchoonschipSettings::default_network())
            .schoonschip()
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        NetworkSchoonschip { settings }.run::<EXPANDSUMS, Aind>(*self)
    }
}
