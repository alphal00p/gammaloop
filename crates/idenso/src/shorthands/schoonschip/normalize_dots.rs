use std::sync::LazyLock;

use spenso::{
    dot, dualizable_, g,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    rank1_, self_dual_,
    structure::{
        abstract_index::AbstractIndex,
        representation::RepName,
        slot::{IsAbstractSlot, SlotMatch, SlotMatcher},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol, representation::FunView},
    function,
    id::Replacement,
};

use crate::{W_, rep_symbols::RS, shorthands::metric::not_slot};

static METRIC_DOT_PRODUCT: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
    let self_dual_stripped1 = self_dual_!(0; W_.d_);
    let self_dual_stripped2 = self_dual_!(0; W_.e_);
    let dualizable_stripped1 = dualizable_!(0; W_.d_);
    let dualizable_stripped2 = dualizable_!(0; W_.e_);

    [
        Replacement::new(
            g!(
                rank1_!(0; W_.d___, &self_dual_stripped1),
                function!(RS.f_, RS.a___, &self_dual_stripped2)
            )
            .to_pattern(),
            dot!(
                rank1_!(0; W_.d___, self_dual_stripped1),
                function!(RS.f_, RS.a___, self_dual_stripped2)
            ),
        )
        .when(not_slot(RS.a___)),
        Replacement::new(
            g!(
                rank1_!(0; W_.d___, &dualizable_stripped1),
                function!(RS.f_, RS.a___, &dualizable_stripped2)
            )
            .to_pattern(),
            dot!(
                rank1_!(0; W_.d___, dualizable_stripped1),
                function!(RS.f_, RS.a___, dualizable_stripped2)
            ),
        )
        .when(not_slot(RS.a___)),
    ]
});

enum DotRewrite {
    Rewritten(Atom),
    Opaque,
    Descend,
}

pub(crate) struct DotNormalizer {
    slots: SlotMatcher,
    rank_one: &'static str,
    metric: Symbol,
    exposes_product: bool,
    metrics: bool,
    rank_one_powers: bool,
    representations: Option<Vec<spenso::structure::representation::LibraryRep>>,
}

/// An odd tensor power exposes a factor which a subsequent contraction
/// can consume. Even powers and nested-vector identities only remove indices.
pub(super) struct DotNormalization {
    pub(super) expression: Atom,
    pub(super) exposes_product: bool,
}

impl DotNormalizer {
    #[cfg(any(test, feature = "reference-cases"))]
    pub(crate) fn run(view: AtomView<'_>) -> Atom {
        Self::normalize(view).expression
    }

    /// Convert compact scalar-product spelling only. Explicit indices and powers
    /// remain available for the structural contractor.
    pub(crate) fn notation(view: AtomView<'_>) -> Atom {
        let compact = Self::with_settings(
            view,
            crate::tensor::ContractSettings {
                metrics: false,
                rank_one: false,
                ..Default::default()
            },
        );
        Self::metric_shorthand_to_dot(compact.as_view())
    }

    pub(crate) fn with_settings(
        view: AtomView<'_>,
        settings: crate::tensor::ContractSettings<'_>,
    ) -> Atom {
        Self::normalize_with_settings(view, settings).expression
    }

    #[cfg(any(test, feature = "reference-cases"))]
    pub(super) fn normalize(view: AtomView<'_>) -> DotNormalization {
        Self::normalize_with_settings(view, Default::default())
    }

    pub(super) fn normalize_with_settings(
        view: AtomView<'_>,
        settings: crate::tensor::ContractSettings<'_>,
    ) -> DotNormalization {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::OutputNormalization,
        );
        let mut normalizer = Self {
            slots: SlotMatcher::default(),
            rank_one: &T.rank1,
            metric: ETS.metric,
            exposes_product: false,
            metrics: settings.metrics,
            rank_one_powers: settings.rank_one,
            representations: settings.representations.map(<[_]>::to_vec),
        };
        let expression = normalizer.apply(view);
        DotNormalization {
            expression,
            exposes_product: normalizer.exposes_product,
        }
    }

    fn apply(&mut self, view: AtomView<'_>) -> Atom {
        match self.normalize_node(view) {
            DotRewrite::Rewritten(result) => return result,
            DotRewrite::Opaque => return view.to_owned(),
            DotRewrite::Descend => {}
        }
        // replace_map treats even an unchanged pruning assignment as a change,
        // rebuilding its ancestors. A visitor can skip opaque payloads without
        // that work when no nested-vector or power identity applies.
        let Some((first_position, first_result)) = self.first_rewrite(view) else {
            return view.to_owned();
        };
        let mut first_result = Some(first_result);
        let mut position = 0;
        view.replace_map(|atom, context, out| {
            let current = position;
            position += 1;
            // Both walks use the same top-down order and opaque-node pruning.
            // The ordinal identifies one occurrence, not every equal subtree.
            if current == first_position {
                **out = first_result.take().unwrap();
            } else if context.parent_type.is_some() {
                match self.normalize_node(atom) {
                    DotRewrite::Rewritten(result) => **out = result,
                    DotRewrite::Opaque => out.set_from_view(&atom),
                    DotRewrite::Descend => {}
                }
            }
        })
    }

    fn first_rewrite(&mut self, view: AtomView<'_>) -> Option<(usize, Atom)> {
        let mut first = None;
        let mut position = 0;
        view.visitor(&mut |atom| {
            if first.is_some() {
                return false;
            }
            let current = position;
            position += 1;
            // apply already tried the root. Reuse that decline in both walks.
            if current == 0 {
                return true;
            }
            match self.normalize_node(atom) {
                DotRewrite::Rewritten(result) => first = Some((current, result)),
                DotRewrite::Opaque => return false,
                DotRewrite::Descend => {}
            }
            first.is_none()
        });
        first
    }

    fn normalize_node(&mut self, atom: AtomView<'_>) -> DotRewrite {
        // Metric construction owns traces and compact-vector substitution.
        // Only nested vectors and powers need this expression walk. Canonical
        // vectors have one final argument carrying their tensor structure;
        // preceding scalar parameters and slot payloads remain opaque.
        match atom {
            AtomView::Fun(function) if function.get_symbol() == T.dot => self
                .normalize_weighted_dot(function)
                .map_or(DotRewrite::Descend, DotRewrite::Rewritten),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                self.normalize_power(base, exponent)
                    .map_or(DotRewrite::Descend, DotRewrite::Rewritten)
            }
            AtomView::Fun(function) if function.get_symbol().has_tag(self.rank_one) => self
                .normalize_nested_vector(function)
                .map_or(DotRewrite::Opaque, DotRewrite::Rewritten),
            AtomView::Fun(function)
                if function.get_symbol().is_scalar()
                    || !matches!(self.slots.classify(atom), SlotMatch::Other) =>
            {
                DotRewrite::Opaque
            }
            _ => DotRewrite::Descend,
        }
    }

    fn normalize_weighted_dot(&mut self, function: FunView<'_>) -> Option<Atom> {
        if !function
            .iter()
            .any(|argument| matches!(argument, AtomView::Mul(_)))
        {
            return None;
        }
        // Reuse the materializer's compatibility and scalar-interface proofs.
        // This only moves coefficients out of an already explicit inner product;
        // it neither pairs indices nor distributes sums inside either operand.
        let (left, right, _, _) = self
            .slots
            .compact_inner_product_parts::<AbstractIndex>(function)?;
        let mut scalar_factors = Vec::new();
        let mut operands = Vec::with_capacity(2);
        for argument in [left, right] {
            if let AtomView::Mul(product) = argument {
                let (position, _) = self
                    .slots
                    .compact_vector_sequence::<AbstractIndex>(product.iter().enumerate(), false)?;
                for (index, factor) in product.iter().enumerate() {
                    if index == position {
                        operands.push(self.apply(factor));
                    } else {
                        scalar_factors.push(self.apply(factor));
                    }
                }
            } else {
                operands.push(self.apply(argument));
            }
        }
        scalar_factors.push(FunctionBuilder::new(T.dot).add_args(operands).finish());
        Some(Atom::mul_many(scalar_factors))
    }

    fn normalize_nested_vector(&mut self, function: FunView<'_>) -> Option<Atom> {
        let inner = self.slots.vector_argument(function)?;
        let AtomView::Fun(inner_function) = inner else {
            return None;
        };
        let compact = self.slots.vector_argument(inner_function)?;
        let representation = self.slots.compact_representation(compact)?;
        if !representation.is_base() {
            return None;
        }
        // The compact final argument identifies even an untagged inner vector;
        // the outer function must have the rank-one tag checked by the caller.
        let outer = Self::replace_last_argument(function, compact.to_owned());
        Some(function!(self.metric, outer, inner))
    }

    fn permits(&mut self, slot: AtomView<'_>) -> bool {
        if self.representations.is_none() {
            return true;
        }
        let Ok(slot) = self.slots.parse::<spenso::structure::representation::LibraryRep, spenso::structure::abstract_index::AbstractIndex>(slot) else {
            return false;
        };
        self.representations.as_ref().is_none_or(|allowed| {
            allowed
                .iter()
                .any(|candidate| candidate.base() == slot.rep().rep.base())
        })
    }

    fn normalize_power(&mut self, base: AtomView<'_>, exponent: AtomView<'_>) -> Option<Atom> {
        let exponent = i64::try_from(exponent).ok()?;
        // Retain the existing integral-power convention: negative even powers
        // normalize, while negative odd and nonintegral powers remain opaque.
        if exponent % 2 == -1 {
            return None;
        }
        let AtomView::Fun(function) = base else {
            return None;
        };
        let square = if self.metrics && function.get_symbol() == self.metric {
            let mut arguments = function.iter();
            if arguments.len() != 2 {
                return None;
            }
            let first_argument = arguments.next().unwrap();
            if !self.permits(first_argument) {
                return None;
            }
            let SlotMatch::Explicit(first) = self.slots.classify(first_argument) else {
                return None;
            };
            let SlotMatch::Explicit(second) = self.slots.classify(arguments.next().unwrap()) else {
                return None;
            };
            if first.is_concrete_index()
                || second.is_concrete_index()
                || !first.representation().is_self_dual()
                || !first.representation().matches(second.representation())
                || first.index() == second.index()
            {
                return None;
            }
            first.dimension().to_owned()
        } else if self.rank_one_powers && function.get_symbol().has_tag(self.rank_one) {
            let argument = self.slots.vector_argument(function)?;
            if !self.permits(argument) {
                return None;
            }
            let SlotMatch::Explicit(slot) = self.slots.classify(argument) else {
                return None;
            };
            if slot.is_concrete_index() || !slot.representation().is_self_dual() {
                return None;
            }
            let compact = Self::replace_last_argument(function, slot.representation().compact());
            function!(self.metric, &compact, &compact)
        } else {
            return None;
        };
        let paired = square.pow(exponent / 2);
        Some(if exponent % 2 == 1 {
            self.exposes_product = true;
            base * paired.as_view()
        } else {
            paired
        })
    }

    fn replace_last_argument(function: FunView<'_>, replacement: Atom) -> Atom {
        FunctionBuilder::new(function.get_symbol())
            .add_args(function.iter().take(function.get_nargs() - 1))
            .add_arg(replacement)
            .finish()
    }

    pub(crate) fn metric_shorthand_to_dot(view: AtomView<'_>) -> Atom {
        view.to_owned().replace_multiple(&*METRIC_DOT_PRODUCT)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::{dualizable_dual_, rep_, structure::abstract_index::AIND_SYMBOLS};
    use symbolica::{function, symbol};
    use symbolica_utils::PatternReplacement;

    #[test]
    fn weighted_dot_notation_extracts_only_scalar_coefficients() {
        crate::test_support::test_initialize();
        let p = T.rank_one_tensor_symbol("weighted_dot_p");
        let q = T.rank_one_tensor_symbol("weighted_dot_q");
        let p = function!(p, spenso::mink!(4));
        let q = function!(q, spenso::mink!(4));
        let scalar = (Atom::var(symbol!("weighted_dot_x")) + Atom::one()).pow(20);
        let inner = dot!(&p, &q);
        let input = dot!(&scalar * &inner * &p, Atom::num(3) * &q);
        let expected = Atom::num(3) * &scalar * inner.pow(2);
        assert_eq!(DotNormalizer::notation(input.as_view()), expected);
        assert_eq!(DotNormalizer::notation(expected.as_view()), expected);

        // An explicit tensor factor is not a scalar coefficient, and notation
        // conversion must not contract it against a repeated index elsewhere.
        let indexed = function!(
            T.rank_one_tensor_symbol("weighted_dot_r"),
            spenso::mink!(4, 1)
        );
        let other = function!(
            T.rank_one_tensor_symbol("weighted_dot_s"),
            spenso::mink!(4, 1)
        );
        let coefficient = &indexed * &other;
        let input = dot!(&coefficient * &p, &q);
        assert_eq!(
            DotNormalizer::notation(input.as_view()),
            &coefficient * dot!(&p, &q)
        );
        let input = dot!(&indexed * &p, &q) * &indexed;
        assert_eq!(DotNormalizer::notation(input.as_view()), input);
    }

    #[test]
    fn first_dot_rewrite_is_constructed_once_per_occurrence() {
        use std::sync::{Arc, Mutex};

        crate::test_support::test_initialize();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let p = spenso::vector_symbol!(
            "dot_first_rewrite_p",
            norm = move |value, _| observed.lock().unwrap().push(value.to_owned())
        );
        let q = T.rank_one_tensor_symbol("dot_first_rewrite_q");
        let scope = symbol!("dot_first_rewrite_scope");
        let scalar = symbol!("dot_first_rewrite_scalar"; Scalar);
        let nested = function!(p, function!(q, spenso::mink!(4)));
        let power = function!(p, spenso::mink!(4, 73401)).pow(3);
        let opaque = function!(scalar, &nested);
        let input = function!(scope, &opaque, &nested, &nested, &power);

        calls.lock().unwrap().clear();
        let first = DotNormalizer::normalize(nested.as_view());
        let second = DotNormalizer::normalize(nested.as_view());
        let third = DotNormalizer::normalize(power.as_view());
        let transcript = calls.lock().unwrap().clone();
        assert!(!transcript.is_empty());
        assert!(third.exposes_product);
        let expected = function!(
            scope,
            &opaque,
            first.expression,
            second.expression,
            third.expression
        );
        calls.lock().unwrap().clear();
        let result = DotNormalizer::normalize(input.as_view());
        assert_eq!(result.expression, expected);
        assert!(result.exposes_product);
        assert_eq!(*calls.lock().unwrap(), transcript);
        calls.lock().unwrap().clear();
        assert_eq!(
            DotNormalizer::run(result.expression.as_view()),
            result.expression
        );
        assert!(calls.lock().unwrap().is_empty());
    }

    #[test]
    fn normalization_reports_new_product_sources() {
        crate::test_support::test_initialize();
        let vector = T.rank_one_tensor_symbol("dot_pending_vector");
        let vector = function!(vector, spenso::mink!(4, 73219));
        let metric = g!(spenso::mink!(4, 73219), spenso::mink!(4, 73223));
        for base in [vector, metric] {
            for exponent in [-4, -3, 2, 3, 4, 5] {
                let input = base.clone().pow(exponent);
                let result = DotNormalizer::normalize(input.as_view());
                assert_eq!(result.expression, DotNormalizer::run(input.as_view()));
                assert_eq!(result.exposes_product, exponent > 1 && exponent % 2 == 1);
                let nested = Atom::var(symbol!("dot_pending_scalar")) + &input;
                assert_eq!(
                    DotNormalizer::normalize(nested.as_view()).exposes_product,
                    result.exposes_product
                );
            }
        }
    }

    static ASYMMETRIC_SCHOONSCHIP_VECTOR_IN_VECTOR: LazyLock<[Replacement; 1]> =
        LazyLock::new(|| {
            let stripped = rep_!(0; W_.d_);
            // A stripped representation identifies the inner vector even when its
            // Symbolica head was not registered with a rank-one tensor tag.
            let inner = function!(RS.f_, RS.a___, &stripped);

            [
                //  p(...,q(..,rep)) is asymmetric, so we replace it with a dot product, using a schoonschiped metric:
                //  p(...,q(..,rep)) => g(p(...,rep), q(..,rep))
                Replacement::new(
                    rank1_!(0; W_.c___, &inner).to_pattern(),
                    g!(rank1_!(0; W_.c___, &stripped), inner,),
                )
                .when(not_slot(RS.a___))
                .min_level(0)
                .max_level(0)
                .level_is_tree_depth(true),
            ]
        });

    static REDUNDANT_METRIC_SCHOONSCHIPS: LazyLock<[Replacement; 4]> = LazyLock::new(|| {
        let self_dual = self_dual_!(0; W_.d_, W_.i_);
        let self_dual_stripped = self_dual_!(0; W_.d_);
        let dualizable = dualizable_!(0; W_.d_, W_.i_);
        let dualizable_stripped = dualizable_!(0; W_.d_);
        let dualizable_dual = dualizable_dual_!(0; W_.d_, W_.i_);

        [
            // g(mu,p(...,rep)) is redundant:
            // g(mu,p(...,rep)) => p(...,mu)
            Replacement::new(
                g!(&self_dual, rank1_!(0; W_.c___, self_dual_stripped)).to_pattern(),
                rank1_!(0; W_.c___, self_dual.clone()),
            ),
            // Same thing but for dualizable
            Replacement::new(
                g!(
                    &dualizable,
                    rank1_!(0; W_.c___, dualizable_stripped.clone())
                )
                .to_pattern(),
                rank1_!(0; W_.c___, dualizable),
            ),
            // Same thing but for dual dualizable
            Replacement::new(
                g!(&dualizable_dual, rank1_!(0; W_.c___, dualizable_stripped)).to_pattern(),
                rank1_!(0; W_.c___, dualizable_dual),
            ),
            // g(mu,p(...)) is also allowed (although a bit weird), here we normalize it to p(...,mu)
            Replacement::new(
                g!(&self_dual, rank1_!(0; W_.c___)).to_pattern(),
                rank1_!(0; W_.c___, self_dual),
            ),
        ]
        .map(|replacement| {
            replacement
                .min_level(0)
                .max_level(0)
                .level_is_tree_depth(true)
        })
    });

    static VECTOR_POWER_NORMALIZATIONS: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
        let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
        let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
        let self_dual_vector = T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual]);
        let self_dual_square = g!(
            T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
            T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
        );

        [
            // Normalize even powers of a vector p(...,mu)^2n -> g(p(...,rep),p(...,rep))^(n/2)
            Replacement::new(
                self_dual_vector.clone().pow(Atom::var(W_.n_)).to_pattern(),
                self_dual_square.pow(Atom::var(W_.n_) / 2),
            )
            .when(W_.n_.filter(even_power)),
            // Normalize odd powers of a vector p(...,mu)^(2n+1) -> p(...,mu) * g(p(...,rep),p(...,rep))^(n/2)
            Replacement::new(
                self_dual_vector.clone().pow(Atom::var(W_.n_)).to_pattern(),
                self_dual_square.pow((Atom::var(W_.n_) - 1) / 2) * self_dual_vector,
            )
            .when(W_.n_.filter(odd_power)),
        ]
    });

    static METRIC_POWER_NORMALIZATIONS: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
        let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
        let self_dual_j = T.self_dual_::<0, _>([W_.d_, W_.j_]);
        let self_dual_metric = g!(&self_dual, &self_dual_j);

        [
            // Normalize even powers of a metric g(rep(dim,i),rep(dim,j))^2n -> dim^(n/2)
            Replacement::new(
                self_dual_metric.pow(Atom::var(W_.n_)).to_pattern(),
                Atom::var(W_.d_).pow(Atom::var(W_.n_) / 2),
            )
            .when(W_.n_.filter(even_power)),
            // Normalize odd powers of a metric g(rep(dim,i),rep(dim,j))^(2n+1) -> dim^(n/2) * g(rep(dim,i),rep(dim,j))
            Replacement::new(
                self_dual_metric.pow(Atom::var(W_.n_)).to_pattern(),
                Atom::var(W_.d_).pow((Atom::var(W_.n_) - 1) / 2) * &self_dual_metric,
            )
            .when(W_.n_.filter(odd_power)),
        ]
    });

    static METRIC_TRACE_NORMALIZATIONS: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
        let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
        let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
        let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

        [
            // g(i,i) -> d
            Replacement::new(g!(&self_dual, &self_dual).to_pattern(), Atom::var(W_.d_)),
            // g(i,dind(i)) -> d
            Replacement::new(
                g!(&dualizable, &dualizable_dual).to_pattern(),
                Atom::var(W_.d_),
            ),
        ]
    });

    fn even_power(exp: AtomView<'_>) -> bool {
        matches!(i64::try_from(exp), Ok(exp) if exp % 2 == 0)
    }

    fn odd_power(exp: AtomView<'_>) -> bool {
        matches!(i64::try_from(exp), Ok(exp) if exp % 2 == 1)
    }

    fn reference(view: AtomView<'_>) -> Atom {
        let asymmetric = ASYMMETRIC_SCHOONSCHIP_VECTOR_IN_VECTOR
            .clone()
            .map(|rule| rule.min_level(0).max_level(None).level_is_tree_depth(false));
        let redundant = REDUNDANT_METRIC_SCHOONSCHIPS
            .clone()
            .map(|rule| rule.min_level(0).max_level(None).level_is_tree_depth(false));
        view.to_owned()
            .replace_multiple(&asymmetric)
            .replace_multiple_repeat(&redundant)
            .replace_multiple(&*VECTOR_POWER_NORMALIZATIONS)
            .replace_multiple(&*METRIC_POWER_NORMALIZATIONS)
            .replace_multiple(&*METRIC_TRACE_NORMALIZATIONS)
    }

    #[test]
    fn direct_dot_normalization_matches_canonical_pattern_oracle() {
        crate::test_support::test_initialize();
        let (d, x, opaque) = symbol!("direct_dot_d", "direct_dot_x", "direct_dot_scope");
        let p = T.rank_one_tensor_symbol("direct_dot_p");
        let q = T.rank_one_tensor_symbol("direct_dot_q");
        // Recognition must not depend on registration in LibraryRep.
        let rep = T.self_dual_symbol("direct_dot_rep");
        let generic_rep = T.representation_symbol("direct_dot_generic_rep");
        let slot = function!(rep, d, 1);
        let other_slot = function!(rep, d, 2);
        let stripped = function!(rep, d);
        let vector = function!(p, x, &slot);
        let metric = g!(&slot, &other_slot);
        let trace = g!(&slot, &slot);
        let nested = function!(p, x, function!(q, 2, &stripped));
        let closed = g!(&slot, function!(p, x, &stripped));
        let mut cases = vec![
            Atom::Zero,
            Atom::var(x),
            Atom::var(p),
            function!(p, &stripped),
            vector.clone(),
            nested.clone(),
            function!(p, function!(opaque, &stripped)),
            function!(p, function!(q, function!(generic_rep, d))),
            closed.clone(),
            metric.clone(),
            trace.clone(),
            g!(function!(p, &stripped), function!(q, &stripped)),
        ];
        for exponent in
            (-4..=4)
                .map(Atom::num)
                .chain([Atom::var(x), Atom::num(1) / 2, Atom::num(-3) / 2])
        {
            for base in [&vector, &metric, &trace, &nested, &closed] {
                cases.push(base.pow(&exponent));
            }
        }
        for (index, case) in cases.iter().enumerate() {
            for expression in [
                case.clone(),
                function!(opaque, case),
                case + Atom::var(x),
                case * (&nested + Atom::var(x)),
                (case + &closed).pow(3),
            ] {
                let normalized = DotNormalizer::run(expression.as_view());
                assert_eq!(
                    normalized,
                    reference(expression.as_view()),
                    "fixture {index}: {expression}"
                );
                assert_eq!(
                    DotNormalizer::run(normalized.as_view()),
                    normalized,
                    "fixture {index}"
                );
            }
        }
        // Traces are scalars before a power is constructed, not metric powers.
        assert_eq!(
            DotNormalizer::run(trace.pow(2).as_view()),
            Atom::var(d).pow(2)
        );
        assert_eq!(
            DotNormalizer::run(trace.pow(3).as_view()),
            Atom::var(d).pow(3)
        );
        assert_eq!(DotNormalizer::run(vector.pow(-3).as_view()), vector.pow(-3));
        let coefficient =
            Atom::add_many((0..256).map(|i| function!(opaque, i)).collect::<Vec<_>>()).pow(7)
                * (Atom::var(x) + function!(opaque, x)).pow(5);
        let expression = &coefficient * (&closed + &nested);
        assert_eq!(
            DotNormalizer::run(expression.as_view()),
            coefficient * DotNormalizer::run((&closed + &nested).as_view())
        );
    }

    #[test]
    fn dot_normalization_requires_canonical_vectors_and_preserves_payloads() {
        crate::test_support::test_initialize();
        let (d, e, x) = symbol!("strict_dot_d", "strict_dot_e", "strict_dot_x");
        let p = T.rank_one_tensor_symbol("strict_dot_p");
        let q = T.rank_one_tensor_symbol("strict_dot_q");
        let scalar = symbol!("strict_dot_scalar"; Scalar);
        let rep = T.self_dual_symbol("strict_dot_rep");
        let other = T.self_dual_symbol("strict_dot_other_rep");
        let dual = T.dualizable_symbol("strict_dot_dual_rep");
        let compact = function!(rep, d);
        let slot = function!(rep, d, 1);
        let nested = function!(q, function!(p, &compact));
        let malformed = [
            function!(p, function!(q, &compact), x),
            function!(p, function!(q, &compact, x)),
            function!(p, function!(q, &compact, &compact)),
            function!(p, function!(q, function!(rep, d, 1, 2))),
            function!(p, function!(scalar, &compact)),
            function!(p, &slot, function!(q, &compact)),
            function!(p, &slot, &slot).pow(2),
            function!(p, &slot, x).pow(2),
            function!(p, function!(dual, d, 1)).pow(2),
            g!(&slot, function!(rep, e, 2)).pow(2),
            g!(&slot, function!(other, d, 2)).pow(2),
            g!(&slot, function!(p, x)),
            g!(&slot, function!(p, function!(rep, e))),
            function!(p, &nested, &slot),
            function!(p, function!(rep, d, &compact)),
            function!(p, function!(rep, &compact)),
            function!(rep, &nested, 1),
            function!(rep, d, &nested),
            function!(AIND_SYMBOLS.dind, &slot, &nested),
            function!(scalar, &nested),
        ];
        // Non-final/multiple slots and missing/mismatched compact slots are no
        // longer repaired by appending an index. Scalar metadata is opaque.
        for expression in malformed {
            assert_eq!(
                DotNormalizer::run(expression.as_view()),
                expression,
                "{expression}"
            );
        }
        // Powers of fixed components do not imply an Einstein sum.
        for index in [AIND_SYMBOLS.cind, AIND_SYMBOLS.find] {
            let first = function!(rep, d, function!(index, 1));
            let second = function!(rep, d, function!(index, 2));
            for exponent in [2, 3, -2] {
                for expression in [
                    function!(p, &first).pow(exponent),
                    g!(&first, &second).pow(exponent),
                    g!(&first, &slot).pow(exponent),
                ] {
                    assert_eq!(DotNormalizer::run(expression.as_view()), expression);
                }
            }
        }
        // Keep exact compound dimensions and arbitrarily large index payloads.
        let dimension = Atom::var(d) + 1;
        let left = function!(rep, &dimension, 1);
        let right = function!(rep, &dimension, 4294967297i64);
        assert_eq!(
            DotNormalizer::run(g!(&left, &right).pow(2).as_view()),
            dimension
        );
        let vector = function!(p, x, &right);
        let stripped = function!(p, x, function!(rep, &dimension));
        assert_eq!(
            DotNormalizer::run(vector.pow(2).as_view()),
            g!(&stripped, &stripped)
        );
    }

    #[test]
    fn dot_normalization_preserves_pruned_ancestor_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicBool, AtomicUsize, Ordering},
        };

        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let enabled = Arc::new(AtomicBool::new(false));
        let observed = Arc::clone(&calls);
        let replace = Arc::clone(&enabled);
        let head = spenso::tensor_symbol!(
            "dot_pruned_ancestor_callback",
            norm = move |_, output| {
                observed.fetch_add(1, Ordering::Relaxed);
                if replace.load(Ordering::Relaxed) {
                    **output = Atom::num(17);
                }
            }
        );
        let p = T.rank_one_tensor_symbol("dot_pruned_callback_p");
        let q = T.rank_one_tensor_symbol("dot_pruned_callback_q");
        let compact = spenso::mink!(4);
        let vector = function!(q, &compact);
        let inert = function!(head, g!(&vector, &vector));
        let work = function!(p, spenso::mink!(4, 1)).pow(2);
        let expression = &inert + &work;
        let normalized = function!(p, &compact);
        let normalized = g!(&normalized, &normalized);

        for replace in [false, true] {
            enabled.store(replace, Ordering::Relaxed);
            calls.store(0, Ordering::Relaxed);
            assert_eq!(DotNormalizer::run(inert.as_view()), inert);
            assert_eq!(calls.load(Ordering::Relaxed), 0);

            let result = DotNormalizer::run(expression.as_view());
            assert_eq!(calls.load(Ordering::Relaxed), 1);
            let unchanged_branch = if replace {
                Atom::num(17)
            } else {
                inert.clone()
            };
            assert_eq!(result, unchanged_branch + &normalized);
        }
    }
}
