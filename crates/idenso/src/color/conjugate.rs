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
            coaf.to_symbolic([RS.a__]).to_pattern(),
            cof.to_symbolic([RS.a__]),
        ),
        Replacement::new(
            cof.to_symbolic([RS.a__]).to_pattern(),
            coaf.to_symbolic([RS.a__]),
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
