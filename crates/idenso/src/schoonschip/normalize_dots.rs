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

use crate::W_;

static ASYMMETRIC_SCHOONSCHIP_VECTOR_IN_VECTOR: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    let stripped = T.rep_::<0, _>([W_.d_]);

    [
        //  p(...,q(..,rep)) is asymmetric, so we replace it with a dot product, using a schoonschiped metric:
        //  p(...,q(..,rep)) => g(p(...,rep), q(..,rep))
        Replacement::new(
            T.rank1_::<0, _>([
                Atom::var(W_.c___),
                T.rank1_::<1, _>([&Atom::var(W_.a___), &stripped]),
            ])
            .to_pattern(),
            ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &stripped]),
                T.rank1_::<1, _>([&Atom::var(W_.a___), &stripped]),
            ),
        ),
    ]
});

static REDUNDANT_METRIC_SCHOONSCHIPS: LazyLock<[Replacement; 4]> = LazyLock::new(|| {
    let index_cond = T.index_fiter(W_.i_);
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

    [
        // g(mu,p(...,rep)) is redundant:
        // g(mu,p(...,rep)) => p(...,mu)
        Replacement::new(
            function!(
                ETS.metric,
                &self_dual,
                T.rank1_::<0, _>([Atom::var(W_.c___), self_dual_stripped])
            )
            .to_pattern(),
            T.rank1_::<0, _>([Atom::var(W_.c___), self_dual.clone()]),
        )
        .with_conditions(index_cond.clone()),
        // Same thing but for dualizable
        Replacement::new(
            function!(
                ETS.metric,
                &dualizable,
                T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_stripped.clone()])
            )
            .to_pattern(),
            T.rank1_::<0, _>([Atom::var(W_.c___), dualizable]),
        )
        .with_conditions(index_cond.clone()),
        // Same thing but for dual dualizable
        Replacement::new(
            function!(
                ETS.metric,
                &dualizable_dual,
                T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_stripped])
            )
            .to_pattern(),
            T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_dual]),
        )
        .with_conditions(index_cond.clone()),
        // g(mu,p(...)) is also allowed (although a bit weird), here we normalize it to p(...,mu)
        Replacement::new(
            function!(
                ETS.metric,
                &self_dual,
                T.rank1_::<0, _>([Atom::var(W_.c___)])
            )
            .to_pattern(),
            T.rank1_::<0, _>([Atom::var(W_.c___), self_dual]),
        )
        .with_conditions(index_cond),
    ]
});

static VECTOR_POWER_NORMALIZATIONS: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
    let index_cond = T.index_fiter(W_.i_);
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
    let self_dual_vector = T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual]);
    let self_dual_square = ETS.metric(
        T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
        T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
    );

    [
        // Normalize even powers of a vector p(...,mu)^2n -> g(p(...,rep),p(...,rep))^(n/2)
        Replacement::new(
            self_dual_vector.clone().pow(Atom::var(W_.n_)).to_pattern(),
            self_dual_square.pow(Atom::var(W_.n_) / 2),
        )
        .with_conditions(index_cond.clone() & W_.n_.filter_single(DotNormalizer::even_power)),
        // Normalize odd powers of a vector p(...,mu)^(2n+1) -> p(...,mu) * g(p(...,rep),p(...,rep))^(n/2)
        Replacement::new(
            self_dual_vector.clone().pow(Atom::var(W_.n_)).to_pattern(),
            self_dual_square.pow((Atom::var(W_.n_) - 1) / 2) * self_dual_vector,
        )
        .with_conditions(index_cond & W_.n_.filter_single(DotNormalizer::odd_power)),
    ]
});

static METRIC_POWER_NORMALIZATIONS: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
    let index_cond = T.index_fiter(W_.i_);
    let other_index_cond = T.index_fiter(W_.j_);
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_j = T.self_dual_::<0, _>([W_.d_, W_.j_]);
    let self_dual_metric = function!(ETS.metric, &self_dual, &self_dual_j);

    [
        // Normalize even powers of a metric g(rep(dim,i),rep(dim,j))^2n -> dim^(n/2)
        Replacement::new(
            self_dual_metric.pow(Atom::var(W_.n_)).to_pattern(),
            Atom::var(W_.d_).pow(Atom::var(W_.n_) / 2),
        )
        .with_conditions(
            index_cond.clone()
                & other_index_cond.clone()
                & W_.n_.filter_single(DotNormalizer::even_power),
        ),
        // Normalize odd powers of a metric g(rep(dim,i),rep(dim,j))^(2n+1) -> dim^(n/2) * g(rep(dim,i),rep(dim,j))
        Replacement::new(
            self_dual_metric.pow(Atom::var(W_.n_)).to_pattern(),
            Atom::var(W_.d_).pow(Atom::var(W_.n_) / 2) * &self_dual_metric,
        )
        .with_conditions(
            index_cond.clone()
                & other_index_cond.clone()
                & W_.n_.filter_single(DotNormalizer::odd_power),
        ),
    ]
});

static METRIC_TRACE_NORMALIZATIONS: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
    let index_cond = T.index_fiter(W_.i_);
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

    [
        // g(i,i) -> d
        Replacement::new(
            function!(ETS.metric, &self_dual, &self_dual).to_pattern(),
            Atom::var(W_.d_),
        )
        .with_conditions(index_cond.clone()),
        // g(i,dind(i)) -> d
        Replacement::new(
            function!(ETS.metric, &dualizable, &dualizable_dual).to_pattern(),
            Atom::var(W_.d_),
        )
        .with_conditions(index_cond),
    ]
});

pub(super) struct DotNormalizer;

impl DotNormalizer {
    pub(super) fn run(view: AtomView<'_>) -> Atom {
        view.to_owned()
            .replace_multiple(&*ASYMMETRIC_SCHOONSCHIP_VECTOR_IN_VECTOR)
            .replace_multiple_repeat(&*REDUNDANT_METRIC_SCHOONSCHIPS)
            .replace_multiple(&*VECTOR_POWER_NORMALIZATIONS)
            .replace_multiple(&*METRIC_POWER_NORMALIZATIONS)
            .replace_multiple(&*METRIC_TRACE_NORMALIZATIONS)
    }

    fn even_power(exp: AtomView<'_>) -> bool {
        matches!(i64::try_from(exp), Ok(exp) if exp % 2 == 0)
    }

    fn odd_power(exp: AtomView<'_>) -> bool {
        matches!(i64::try_from(exp), Ok(exp) if exp % 2 == 1)
    }
}
