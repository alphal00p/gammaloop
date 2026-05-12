use std::sync::LazyLock;

use spenso::{
    s,
    structure::representation::{Minkowski, RepName, Representation},
};

use crate::representations::{Bispinor, ColorAdjoint, ColorFundamental};

#[derive(Clone)]
pub(crate) struct TestReps {
    pub mink4: Representation<Minkowski>,
    pub mink_d: Representation<Minkowski>,
    pub bis4: Representation<Bispinor>,
    pub bis_d: Representation<Bispinor>,
    pub cof_nc: Representation<ColorFundamental>,
    pub coad_na: Representation<ColorAdjoint>,
}

pub(crate) static TEST_REPS: LazyLock<TestReps> = LazyLock::new(|| {
    TestReps::initialize_symbols();
    TestReps::build()
});

pub(crate) fn test_initialize() -> &'static TestReps {
    &TEST_REPS
}

impl TestReps {
    pub(crate) fn new() -> Self {
        test_initialize().clone()
    }

    fn build() -> Self {
        Self {
            mink4: Minkowski {}.new_rep(4),
            mink_d: Minkowski {}.new_rep(s!(d)),
            bis4: Bispinor {}.new_rep(4),
            bis_d: Bispinor {}.new_rep(s!(d)),
            cof_nc: ColorFundamental {}.new_rep(s!(Nc)),
            coad_na: ColorAdjoint {}.new_rep(s!(NA)),
        }
    }

    fn initialize_symbols() {
        crate::representations::initialize();

        let tags = &spenso::network::tags::SPENSO_TAG;
        let _ = tags.rank_one_tensor_symbol("P");
        let _ = tags.rank_one_tensor_symbol("Q");
        let _ = tags.rank_one_tensor_symbol("K");

        let _ = spenso::p!();
        let _ = spenso::vector_symbol!(q);
        let _ = spenso::vector_symbol!(P);
        let _ = spenso::vector_symbol!(Q);
        let _ = spenso::vector_symbol!(K);
    }
}
