use std::sync::LazyLock;

use spenso::structure::representation::RepName;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    id::Replacement,
};

use crate::{color_t, rep_symbols::RS, representations::ColorFundamental};

static COLOR_CONJ_GENERATOR_TRANSPOSITIONS: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    let cof = ColorFundamental {};
    let coaf = ColorFundamental {}.dual();

    [Replacement::new(
        color_t!(
            RS.i__,
            cof.to_symbolic([RS.d_, RS.a_]),
            coaf.to_symbolic([RS.d_, RS.b_]),
        )
        .to_pattern(),
        color_t!(
            RS.i__,
            coaf.to_symbolic([RS.d_, RS.b_]),
            cof.to_symbolic([RS.d_, RS.a_]),
        ),
    )]
});

static COLOR_CONJ_REPRESENTATION_SWAPS: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
    let cof = ColorFundamental {};
    let coaf = ColorFundamental {}.dual();

    [
        Replacement::new(
            coaf.to_symbolic([RS.d_, RS.a_]).to_pattern(),
            cof.to_symbolic([RS.d_, RS.a_]),
        ),
        Replacement::new(
            cof.to_symbolic([RS.d_, RS.a_]).to_pattern(),
            coaf.to_symbolic([RS.d_, RS.a_]),
        ),
    ]
});

struct ColorConjugator;

impl ColorConjugator {
    fn run(expression: AtomView<'_>) -> Atom {
        expression
            .to_owned()
            .replace_multiple(&*COLOR_CONJ_GENERATOR_TRANSPOSITIONS)
            .replace_multiple(&*COLOR_CONJ_REPRESENTATION_SWAPS)
    }
}

pub fn color_conj_impl(expression: AtomView<'_>) -> Atom {
    ColorConjugator::run(expression)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        IndexTooling, color::ColorSimplifier, shorthands::metric::MetricSimplifier,
        test_support::test_initialize,
    };
    use spenso::structure::abstract_index::AbstractIndex;
    use symbolica::parse;

    #[test]
    fn adjoint_color_generator_has_positive_su3_norm() {
        test_initialize();
        let generator = parse!(
            "t(coad(8,a),cof(3,i),dind(cof(3,j)))",
            default_namespace = "spenso"
        );
        let adjoint = generator.dirac_adjoint::<AbstractIndex>().unwrap();
        assert_eq!(
            adjoint,
            parse!(
                "t(coad(8,a),cof(3,j),dind(cof(3,i)))",
                default_namespace = "spenso"
            )
        );
        let norm = (generator * adjoint)
            .simplify_color()
            .to_cof_dimension_invariants()
            .simplify_metrics();
        assert_eq!(norm, Atom::num(4));
    }

    #[test]
    fn color_conjugation_dualizes_open_network_and_is_involutive() {
        test_initialize();
        let network = Atom::i()
            * parse!(
                "t(coad(8,a),cof(3,i),dind(cof(3,j)))*g(cof(3,j),dind(cof(3,k)))",
                default_namespace = "spenso"
            );
        let conjugate = network.spenso_conj();
        let expected = -Atom::i()
            * parse!(
                "t(coad(8,a),cof(3,j),dind(cof(3,i)))*g(cof(3,k),dind(cof(3,j)))",
                default_namespace = "spenso"
            );
        assert_eq!(conjugate, expected);
        assert_eq!(conjugate.spenso_conj(), network);
    }

    #[test]
    fn color_conjugation_preserves_scalar_representation_invariants() {
        test_initialize();
        let invariant = parse!("cas(2,cof(3))", default_namespace = "spenso");
        assert_eq!(invariant.spenso_conj(), invariant);
    }
}
