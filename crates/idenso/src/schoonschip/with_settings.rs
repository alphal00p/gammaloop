use std::sync::LazyLock;

use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    tensors::parametric::atomcore::PatternReplacement,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::Replacement,
};

use crate::{W_, metric::not_slot};

use super::{
    api::Schoonschip, chain_like::simplify_chain_like_metric_products,
    settings::SchoonschipSettings,
};

static METRIC_FUNCTION_CONTRACTIONS: LazyLock<[Replacement; 3]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);
    let function_with_replacement = function!(W_.a_, W_.a___, W_.c_, W_.b___);

    [
        Replacement::new(
            (function!(ETS.metric, &self_dual, W_.c_)
                * function!(W_.a_, W_.a___, &self_dual, W_.b___))
            .to_pattern(),
            function_with_replacement.clone(),
        ),
        Replacement::new(
            (function!(ETS.metric, W_.c_, &dualizable)
                * function!(W_.a_, W_.a___, &dualizable_dual, W_.b___))
            .to_pattern(),
            function_with_replacement.clone(),
        ),
        Replacement::new(
            (function!(ETS.metric, W_.c_, &dualizable_dual)
                * function!(W_.a_, W_.a___, &dualizable, W_.b___))
            .to_pattern(),
            function_with_replacement.clone(),
        ),
    ]
});

static DOT_PRODUCT: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);

    [Replacement::new(
        (function!(W_.f_, W_.a___, &self_dual) * function!(W_.g_, W_.b___, &self_dual))
            .to_pattern(),
        ETS.metric(
            function!(W_.f_, W_.a___, &self_dual_stripped),
            function!(W_.g_, W_.b___, &self_dual_stripped),
        ),
    )
    .with_conditions(not_slot(W_.a___) & not_slot(W_.b___))]
});

static SELF_DUAL_POWER: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);

    [Replacement::new(
        function!(W_.f_, W_.a___, &self_dual)
            .pow(Atom::num(2))
            .to_pattern(),
        ETS.metric(
            function!(W_.f_, W_.a___, &self_dual_stripped),
            function!(W_.f_, W_.a___, &self_dual_stripped),
        ),
    )
    .with_conditions(not_slot(W_.a___))]
});

static RANK_ONE_PRODUCT_REPLACEMENTS: LazyLock<[Replacement; 3]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

    [
        Replacement::new(
            (T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual])
                * T.rank1_::<1, _>([&Atom::var(W_.a___), &self_dual]))
            .to_pattern(),
            ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                T.rank1_::<1, _>([&Atom::var(W_.a___), &self_dual_stripped]),
            ),
        ),
        Replacement::new(
            T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual])
                .pow(Atom::num(2))
                .to_pattern(),
            ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
            ),
        ),
        Replacement::new(
            (T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable])
                * T.rank1_::<1, _>([&Atom::var(W_.a___), &dualizable_dual]))
            .to_pattern(),
            ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                T.rank1_::<1, _>([&Atom::var(W_.a___), &dualizable_stripped]),
            ),
        ),
    ]
});

static BARE_PRODUCT_CONTRACTIONS: LazyLock<[Replacement; 3]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

    [
        Replacement::new(
            (function!(W_.a_, W_.a___, &self_dual, W_.b___)
                * T.rank1_::<0, _>([Atom::var(W_.c___), self_dual.clone()]))
            .to_pattern(),
            function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                W_.b___
            ),
        ),
        Replacement::new(
            (function!(W_.a_, W_.a___, &dualizable, W_.b___)
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_dual]))
            .to_pattern(),
            function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                W_.b___
            ),
        ),
        Replacement::new(
            (function!(W_.a_, W_.a___, &dualizable_dual, W_.b___)
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable]))
            .to_pattern(),
            function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_stripped]),
                W_.b___
            ),
        ),
    ]
});

struct SchoonschipWithSettings<'a> {
    settings: &'a SchoonschipSettings,
}

impl SchoonschipWithSettings<'_> {
    fn run(&self, view: AtomView<'_>) -> Atom {
        let mut current = view.to_owned();
        loop {
            let next = self.apply_once(current.as_view());
            if next == current {
                return next;
            }
            current = next;
        }
    }

    fn apply_once(&self, view: AtomView<'_>) -> Atom {
        let metric_simplified = view
            .normalize_dots()
            .to_owned()
            .replace_multiple_repeat(&*METRIC_FUNCTION_CONTRACTIONS);

        let simplified = metric_simplified
            .replace_multiple(&*DOT_PRODUCT)
            .replace_multiple(&*SELF_DUAL_POWER)
            .replace_multiple(&*RANK_ONE_PRODUCT_REPLACEMENTS)
            // Bare product contraction: vector(rep(d,i)) * T(..,rep(d,i),..)
            // becomes T(..,vector(rep(d)),..). Network contraction has a more
            // structured version of this rule and does not use this pass.
            .replace_multiple_repeat(&*BARE_PRODUCT_CONTRACTIONS);

        let simplified = if self.settings.simplify_chain_like_functions {
            simplify_chain_like_metric_products(simplified.as_view())
        } else {
            simplified
        };

        simplified.normalize_dots()
    }
}

pub(super) fn schoonschip_with_settings(
    view: AtomView<'_>,
    settings: &SchoonschipSettings,
) -> Atom {
    SchoonschipWithSettings { settings }.run(view)
}
