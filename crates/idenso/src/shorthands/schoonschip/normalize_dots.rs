use std::sync::LazyLock;

use spenso::{
    dot, dualizable_, dualizable_dual_, g,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    rank1_, rep_, self_dual_,
    structure::abstract_index::AIND_SYMBOLS,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    id::Replacement,
};

use crate::W_;

static ASYMMETRIC_SCHOONSCHIP_VECTOR_IN_VECTOR: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    let stripped = rep_!(0; W_.d_);

    [
        //  p(...,q(..,rep)) is asymmetric, so we replace it with a dot product, using a schoonschiped metric:
        //  p(...,q(..,rep)) => g(p(...,rep), q(..,rep))
        Replacement::new(
            rank1_!(0; W_.c___, rank1_!(1; W_.a___, &stripped)).to_pattern(),
            g!(
                rank1_!(0; W_.c___, &stripped),
                rank1_!(1; W_.a___, &stripped),
            ),
        )
        .level_range((0, Some(0)))
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
            .level_range((0, Some(0)))
            .level_is_tree_depth(true)
    })
});

static METRIC_DOT_PRODUCT: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
    let self_dual_stripped1 = self_dual_!(0; W_.d_);
    let self_dual_stripped2 = self_dual_!(0; W_.e_);
    let dualizable_stripped1 = dualizable_!(0; W_.d_);
    let dualizable_stripped2 = dualizable_!(0; W_.e_);

    [
        Replacement::new(
            g!(
                rank1_!(0; W_.d___, &self_dual_stripped1),
                rank1_!(1; W_.c___, &self_dual_stripped2)
            )
            .to_pattern(),
            dot!(
                rank1_!(0; W_.d___, self_dual_stripped1),
                rank1_!(1; W_.c___, self_dual_stripped2)
            ),
        ),
        Replacement::new(
            g!(
                rank1_!(0; W_.d___, &dualizable_stripped1),
                rank1_!(1; W_.c___, &dualizable_stripped2)
            )
            .to_pattern(),
            dot!(
                rank1_!(0; W_.d___, dualizable_stripped1),
                rank1_!(1; W_.c___, dualizable_stripped2)
            ),
        ),
    ]
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
        .when(W_.n_.filter(DotNormalizer::even_power)),
        // Normalize odd powers of a vector p(...,mu)^(2n+1) -> p(...,mu) * g(p(...,rep),p(...,rep))^(n/2)
        Replacement::new(
            self_dual_vector.clone().pow(Atom::var(W_.n_)).to_pattern(),
            self_dual_square.pow((Atom::var(W_.n_) - 1) / 2) * self_dual_vector,
        )
        .when(W_.n_.filter(DotNormalizer::odd_power)),
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
        .when(W_.n_.filter(DotNormalizer::even_power)),
        // Normalize odd powers of a metric g(rep(dim,i),rep(dim,j))^(2n+1) -> dim^(n/2) * g(rep(dim,i),rep(dim,j))
        Replacement::new(
            self_dual_metric.pow(Atom::var(W_.n_)).to_pattern(),
            Atom::var(W_.d_).pow((Atom::var(W_.n_) - 1) / 2) * &self_dual_metric,
        )
        .when(W_.n_.filter(DotNormalizer::odd_power)),
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

pub(crate) struct DotNormalizer;

impl DotNormalizer {
    pub(super) fn run(view: AtomView<'_>) -> Atom {
        // Traverse once per stage and only invoke the matcher at eligible roots.
        // Ordinary momentum components and scalar products need no wildcard
        // matching. Keep the original top-down traversal and stage ordering:
        // metric powers must be normalized before their traces, for example.
        // Inspect all arguments because symmetric heads can permute their slots.
        let mut current = view.replace_map(|atom, _context, out| {
            let AtomView::Fun(fun) = atom else { return };
            if fun.get_symbol().has_tag(&T.rank1)
                && fun.iter().any(|arg| {
                    matches!(arg, AtomView::Fun(inner)
                        if inner.get_symbol().has_tag(&T.rank1)
                            && inner.iter().any(|arg| matches!(arg, AtomView::Fun(rep)
                                if rep.get_symbol().has_tag(&T.representation))))
                })
            {
                let mut replaced = Atom::new();
                if atom
                    .replace_multiple_into(&*ASYMMETRIC_SCHOONSCHIP_VECTOR_IN_VECTOR, &mut replaced)
                {
                    **out = replaced;
                }
            }
        });
        current.repeat_map(|view| {
            view.replace_map(|atom, _context, out| {
                let AtomView::Fun(fun) = atom else { return };
                if fun.get_symbol() == ETS.metric
                    && fun.iter().any(|arg| {
                        matches!(arg, AtomView::Fun(rep)
                        if rep.get_symbol().has_tag(&T.representation)
                            || rep.get_symbol() == AIND_SYMBOLS.dind)
                    })
                {
                    let mut replaced = Atom::new();
                    if atom.replace_multiple_into(&*REDUNDANT_METRIC_SCHOONSCHIPS, &mut replaced) {
                        **out = replaced;
                    }
                }
            })
        });
        current
            .replace_multiple(&*VECTOR_POWER_NORMALIZATIONS)
            .replace_multiple(&*METRIC_POWER_NORMALIZATIONS)
            .replace_multiple(&*METRIC_TRACE_NORMALIZATIONS)
    }

    pub(crate) fn metric_shorthand_to_dot(view: AtomView<'_>) -> Atom {
        view.to_owned().replace_multiple(&*METRIC_DOT_PRODUCT)
    }

    fn even_power(exp: AtomView<'_>) -> bool {
        matches!(i64::try_from(exp), Ok(exp) if exp % 2 == 0)
    }

    fn odd_power(exp: AtomView<'_>) -> bool {
        matches!(i64::try_from(exp), Ok(exp) if exp % 2 == 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::{function, symbol};
    use symbolica_utils::PatternReplacement;

    #[test]
    fn staged_dot_normalization_matches_unrestricted_pattern_oracle() {
        crate::test_support::test_initialize();
        let asymmetric = ASYMMETRIC_SCHOONSCHIP_VECTOR_IN_VECTOR
            .clone()
            .map(|rule| rule.level_range((0, None)).level_is_tree_depth(false));
        let redundant = REDUNDANT_METRIC_SCHOONSCHIPS
            .clone()
            .map(|rule| rule.level_range((0, None)).level_is_tree_depth(false));
        // Keep the former whole-expression matcher as an independent oracle.
        // Its five stages, including the redundant-metric fixed point, must
        // agree exactly without expanding these diagnostic-only expressions.
        let reference = |view: AtomView<'_>| {
            view.to_owned()
                .replace_multiple(&asymmetric)
                .replace_multiple_repeat(&redundant)
                .replace_multiple(&*VECTOR_POWER_NORMALIZATIONS)
                .replace_multiple(&*METRIC_POWER_NORMALIZATIONS)
                .replace_multiple(&*METRIC_TRACE_NORMALIZATIONS)
        };

        let (d, e, x, opaque) = symbol!(
            "dot_guard_d",
            "dot_guard_e",
            "dot_guard_x",
            "dot_guard_opaque"
        );
        let p = T.rank_one_tensor_symbol("dot_guard_p");
        let q = T.rank_one_tensor_symbol("dot_guard_q");
        let symmetric = symbol!("dot_guard_symmetric"; Symmetric;
            tags = [&T.tensor, &T.rank1]);
        let linear = symbol!("dot_guard_linear"; Linear;
            tags = [&T.tensor, &T.rank1]);
        // These tagged representations are deliberately absent from LibraryRep.
        let rep = T.self_dual_symbol("dot_guard_rep");
        let dualizable = T.dualizable_symbol("dot_guard_dualizable");
        let generic_rep = T.representation_symbol("dot_guard_generic_rep");
        let slot = function!(rep, d, 1);
        let other_slot = function!(rep, d, 2);
        let stripped = function!(rep, d);
        let dual_slot = function!(dualizable, d, 1);
        let dual_stripped = function!(dualizable, d);
        let vector = function!(p, x, &slot);
        let metric = g!(&slot, &other_slot);
        let traced_metric = g!(&slot, &slot);
        let nested = function!(p, x, function!(q, 2, &stripped));
        let closed = g!(&slot, function!(p, x, &stripped));

        let mut cases = vec![
            Atom::num(0),
            Atom::var(x),
            Atom::var(p),
            function!(p, &stripped),
            vector.clone(),
            nested.clone(),
            function!(p, function!(q, function!(p, &stripped))),
            function!(p, function!(opaque, &nested), &slot),
            function!(p, function!(q, function!(generic_rep, d))),
            function!(symmetric, x, function!(q, &stripped), 7),
            function!(p, function!(symmetric, x, &stripped, 7)),
            function!(linear, Atom::num(2) * function!(q, &stripped)),
            function!(p, function!(linear, &stripped)),
            closed.clone(),
            g!(&slot, Atom::var(p)),
            g!(&slot, function!(p, x)),
            g!(&slot, function!(p, x, function!(rep, e))),
            g!(&slot, function!(p, x, function!(generic_rep, d))),
            g!(&slot, function!(opaque, &stripped)),
            g!(&dual_slot, function!(q, &dual_stripped)),
            g!(AIND_SYMBOLS.dual(&dual_slot), function!(q, &dual_stripped)),
            g!(&dual_slot, function!(q, function!(dualizable, e))),
            g!(&dual_slot, AIND_SYMBOLS.dual(&dual_slot)),
            g!(function!(p, &stripped), function!(q, &stripped)),
            g!(&slot, function!(p, &nested)),
            g!(&slot, function!(p, &closed)),
            metric.clone(),
            traced_metric.clone(),
            function!(ETS.metric, &slot),
            function!(ETS.metric, &slot, &other_slot, &nested),
            function!(AIND_SYMBOLS.dind, &slot, &nested),
        ];
        let exponents =
            (-4..=4)
                .map(Atom::num)
                .chain([Atom::var(x), Atom::num(1) / 2, Atom::num(-3) / 2]);
        for exponent in exponents {
            for base in [&vector, &metric, &traced_metric, &nested, &closed] {
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
                assert_eq!(
                    DotNormalizer::run(expression.as_view()),
                    reference(expression.as_view()),
                    "normalization fixture {index}"
                );
            }
        }

        // Trace and negative-power behavior depend on the original stage order
        // and exponent filters, not on a bottom-up algebraic interpretation.
        assert_eq!(
            DotNormalizer::run(traced_metric.pow(2).as_view()),
            Atom::var(d)
        );
        assert_eq!(DotNormalizer::run(vector.pow(-3).as_view()), vector.pow(-3));

        let coefficient =
            Atom::add_many((0..256).map(|i| function!(opaque, i)).collect::<Vec<_>>()).pow(7)
                * (Atom::var(x) + function!(opaque, x)).pow(5);
        let expression = &coefficient * (&closed + &nested);
        let normalized = DotNormalizer::run(expression.as_view());
        assert_eq!(normalized, reference(expression.as_view()));
        assert_eq!(
            normalized,
            coefficient * DotNormalizer::run((&closed + &nested).as_view())
        );
    }
}
